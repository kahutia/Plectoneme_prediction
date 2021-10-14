#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function PlectonemeCode(Swave, Twave, Wwave, Dwave, CwaveRR, CwaveTT)

wave/T Swave		//the sequence of DNA we are considering
wave Twave, Wwave, Dwave, CwaveRR, CwaveTT		//Twist, Wedge, and Direction parameters; Covariance parameters for flexibility (roll-roll and tilt-tilt)

variable ii, Letter1, Letter2, index, alpha_n, beta_n, omDiv2_n

variable SeqLength=dimsize(Swave,0)
variable rise=0.339

make/d/o/n=(SeqLength,4) DNApath, DNApathMajorGroove  //These 2 paths will trace the center and the major groove
make/d/o/n=(SeqLength,2, 2) BasepairCovariance, LocalCovariance  //covariance matrices to estimate local stoffness along the sequence
make/d/o/n=(2, 2) BendRot, Covar, LocalCov  //more matrices for stiffness calculation
DNApath=0
make/d/o/n=(SeqLength) CurvatureSequence
make/d/o/n=(SeqLength) Sequence_phase, Sequence_angle_energy, Sequence_angle_exp, EndEffects  //Sequence_angle_exp will eventually store the weight assigned to a plectoneme at each position on the DNA
CurvatureSequence=0;
Sequence_phase=0

make/d/o/n=(4,4) T_n, Romega_n, Q_n, Rzplus, Rx, Rzminus, Minverse_tot, M_tot
make/d/o/n=4 StartPos, StartPosMG
StartPos={0,0,0,1}
StartPosMG={1,0,0,1}
make/d/o/n=4 OldPos

DNApath[0][]= StartPos[q]

OldPos=StartPos

Minverse_tot={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}
M_tot=Minverse_tot
T_n={{1,0,0,0},{0,1,0,0},{0,0,1,-rise/2},{0,0,0,1}}

//This loop finds the 3D path of the relaxed DNA

Letter1=0;
For (ii=1; ii<SeqLength; ii+=1)
	If  ((char2num(Swave[ii]))==char2num("A")) 
		Letter2=0
	Elseif ((char2num(Swave[ii]))==char2num("C"))
		Letter2=1
	Elseif ((char2num(Swave[ii]))==char2num("G"))
		Letter2=2
	Else
		Letter2=3
	EndIf
	index = 4*Letter1+Letter2			//The index defines the current dinucleotide, AA=0, AC=1,...TT=15
	
	Sequence_phase[ii]=Sequence_phase[ii-1]+Twave[index]					//This is used to measure how far around the DNA the major groove has rotated relative to the first base pair
	BendRot={{cos(Sequence_phase[ii]), sin(Sequence_phase[ii])},{-sin(Sequence_phase[ii]),cos(Sequence_phase[ii])}}  //Rotation matrix 
	Covar={{CwaveRR[index], 0},{0,CwaveTT[index]}}		//The covariance matrix for the current basepair, expressed in the coordinates of the current basepair
	MatrixOp/o CovRot = ( BendRot x Covar x BendRot^t)	//Rotating the covariance matrix so it will line up with its neighbors
	BasepairCovariance[ii][][]=CovRot[q][r]					//Rotated covariance matrix at position is recorded

	//these next steps are outlined in the paper I sent
	omDiv2_n=Twave[index]/2
	
	Romega_n={{cos(omDiv2_n),sin(omDiv2_n),0,0},{-sin(omDiv2_n),cos(omDiv2_n),0,0},{0,0,1,0},{0,0,0,1}}
	alpha_n=Wwave[index]

	beta_n=Dwave[index]-pi/2
	
	Rzplus={{cos(beta_n),sin(beta_n),0,0},{-sin(beta_n),cos(beta_n),0,0},{0,0,1,0},{0,0,0,1}}
	Rx={{1,0,0,0},{0,cos(-alpha_n),sin(-alpha_n),0},{0,-sin(-alpha_n),cos(-alpha_n),0},{0,0,0,1}}
	Rzminus={{cos(-beta_n),sin(-beta_n),0,0},{-sin(-beta_n),cos(-beta_n),0,0},{0,0,1,0},{0,0,0,1}}
	
	MatrixOp/o Q_n = ( Rzminus x Rx x Rzplus )
	MatrixOp/o Minverse_n = Inv( T_n x Romega_n x Q_n x Romega_n x T_n)
	MatrixOp/o Minverse_new = (  Minverse_n x Minverse_tot)

	Minverse_tot=Minverse_new //Updating the total tranformation matrix
	
	MatrixOp/o CurrentPos = ( Minverse_tot^t x StartPos )   //Calculate the coordinates of the current basepair

	MatrixOp/o CurrentPosMG = ( Minverse_tot^t x StartPosMG )
	
	DNApath[ii][]= CurrentPos[q]
	DNApathMajorGroove[ii][] = CurrentPosMG[q]
	
	Letter1=Letter2
Endfor

print "path calculated"

//Make the curvature calculation

variable CurveWindow  //must be even
Sequence_angle_energy=0
Sequence_angle_exp=0
Variable CircFrac=0.667	//We assume the plectoneme tip makes a 240� arc before joining the bulk plectoneme region
Variable BindLength=450   //experimentally, ~450 nt are bound to the surface at each end of the DNA
Variable AvePlecLength=1000	
Variable EnergyOffset
Variable Z1, Z2, Z3, Z4, Z5, Z6, Z7 	//These are used to assign Boltzmann weights to bending in 8 possible directions
Variable Cnorm, E_base

Variable tanlength=10 //must be even. This is number of basepairs used to calculate the local tangent vectors.
Variable VectorMag, CosCurve, SinCurve


For (CurveWindow=40; CurveWindow<120; CurveWindow+=8)
	print CurveWindow
	LocalCovariance=BasepairCovariance
	smooth/B/DIM=0 (CurveWindow+1), LocalCovariance		//Find covariance matrix over the curvature window
	
	make/d/o/n=(SeqLength,3) TanVector, NormVector, CurveVector
	make/d/o/n=(SeqLength) CurveMag, CurvePhase, HalfCurveMag, HalfCurvePhase
	TanVector=0
	NormVector=DNApathMajorGroove-DNApath  //identifies the normal vector alligned with the major groove
	CurveVector=0
	make/d/o/n=3 CurrentTan, CurrentCurve, tP, tM, CurveCross, CurrentNorm  //tangent vectors
		
	// find the tan vectors over tanlength
	For (ii=tanlength/2; ii<SeqLength-tanlength/2; ii+=1)
		CurrentTan[]=DNApath[ii+tanlength/2][p]-DNApath[ii-tanlength/2][p]
		VectorMag=sqrt(MatrixDot(CurrentTan, CurrentTan))
		TanVector[ii][]=CurrentTan[q]/VectorMag  //Normalizes tangent vector to unit length
	EndFor

	// find the curvature vectors and values over curvewindow
	For (ii=CurveWindow/2; ii<SeqLength-CurveWindow/2; ii+=1)
		tP[]=TanVector[ii+CurveWindow/2][p]					//plus tan vector
		tM[]=TanVector[ii-CurveWindow/2][p]					//minus tan vector
		CurrentCurve ={tP[1]*tM[2]-tP[2]*tM[1],tP[2]*tM[0]-tP[0]*tM[2],tP[0]*tM[1]-tP[1]*tM[0]}   //Cross product
		CurveVector[ii][]=CurrentCurve[q]						//curvature vector is recorded at this position
	//	CurveMag[ii]=acos(MatrixDot(tP, tM))
		CurveMag[ii]=asin(sqrt(MatrixDot(CurrentCurve, CurrentCurve)))
		CurrentCurve/=CurveMag[ii]							//normalize the curvature vector to track direction
		
		//Calculates the phase angle of the curvature relative to major groove at start of DNA 	
		CurrentTan[]=TanVector[ii][p]
		CurrentNorm[]=NormVector[ii][p]
		CosCurve=MatrixDot(CurrentCurve, CurrentNorm)
		tP=CurrentNorm	//Only using these variables as placeholders b/c of short names, otherwise cross product calculation makes a long line of code
		tM=CurrentCurve
		CurveCross={tP[1]*tM[2]-tP[2]*tM[1],tP[2]*tM[0]-tP[0]*tM[2],tP[0]*tM[1]-tP[1]*tM[0]} 
		SinCurve=MatrixDot(CurveCross, CurrentTan)
		CurvePhase[ii]=mod(atan2(SinCurve,CosCurve)+Sequence_phase[ii],2*pi)  //curvature phase is recorded.
	EndFor
	
	EnergyOffset=25-CurveWindow*0.334*3/4.06   //adds energy penalty from pulling in DNA ends against a force
	
	For (ii=BindLength; ii<SeqLength-BindLength; ii+=1)
		Covar=LocalCovariance[ii][p][q]
		BendRot={{cos(CurvePhase[ii]), sin(CurvePhase[ii])},{-sin(CurvePhase[ii]),cos(CurvePhase[ii])}}
		MatrixOp/o CovRot = ( BendRot x Covar x BendRot^t)		//local covariance matrix alligned to major groove
		BendRot={{cos(pi/4), sin(pi/4)},{-sin(pi/4),cos(pi/4)}}
		MatrixOp/o CovRot45 = ( BendRot x Covar x BendRot^t)		//local covariance matrix alligned 45� to major groove
		Cnorm=CurveMag[ii]/(2*pi*CircFrac)
//				Cnorm=0				//Uncomment to compare to straight DNA with variable stiffness
		E_base=CircFrac^2*3000/Curvewindow
		Z1=exp(-E_base/CovRot[0][0]*((1-Cnorm)^2)+EnergyOffset)			//Bend in direction of curve
		Z2=exp(-E_base/CovRot[0][0]*((1+Cnorm)^2)+EnergyOffset)			//Bend against curve
		Z3=exp(-E_base/CovRot[1][1]*(1-(Cnorm)^2)+EnergyOffset)			//Bend perpendicular to curve
		Z4=exp(-E_base/CovRot45[0][0]*(sqrt(Cnorm^2/2+1)-Cnorm/sqrt(2))^2+EnergyOffset)		//Bend at 45�, 135�, 225�, and 315�
		Z5=exp(-E_base/CovRot45[0][0]*(sqrt(Cnorm^2/2+1)+Cnorm/sqrt(2))^2+EnergyOffset)
		Z6=exp(-E_base/CovRot45[1][1]*(sqrt(Cnorm^2/2+1)-Cnorm/sqrt(2))^2+EnergyOffset)
		Z7=exp(-E_base/CovRot45[1][1]*(sqrt(Cnorm^2/2+1)+Cnorm/sqrt(2))^2+EnergyOffset)
	
		Sequence_angle_exp[ii]+=Z1+Z2+2*Z3+Z4+Z5+Z6+Z7
	endfor

endfor

Sequence_angle_energy=-ln(Sequence_angle_exp)  //Backing out the implied energy landscape from the summed Boltzmann weights

EndEffects=max(0,min(1,(p-BindLength)/AvePlecLength)*min(1,(SeqLength-p-BindLength)/AvePlecLength))  //takes into account the effects of the handles, including limited plectoneme growth near the attachment points
Sequence_angle_exp*=EndEffects

Duplicate/O Sequence_angle_exp,Sequence_angle_exp_smth
Smooth/E=2/F/B=64 300, Sequence_angle_exp_smth  //this command applies an approximate Gaussian smooth, though in reality it is 64 sequential boxcar smoothing operations
wavestats/q Sequence_angle_exp_smth
Sequence_angle_exp_smth/=V_avg

End
