(* ::Package:: *)

(* Setup of the basic quantities *)
ClearAll["Global`*"]
beginTime=SessionTime[];
notebookDir="/mnt/raid/cieslik/SymV05/b_1";
\[Xi]h=2; (* Location of the horizon *) 
\[Xi]f=3; (* Location of the photon orbit *) 
(* Mass of the test particle and the black hole *)
M = 1; 
m = 1;
(* \[Alpha] and \[Beta] from Maxwell-Juttner distribytion function *)
\[Alpha] = 1;
\[Beta] = 1; 

x\[Xi]0=1000; (* Outer radius of the numerical grid *)

v=0.5; (* Velocity of the boost *)
\[Gamma]=1/Sqrt[1-v^2]; (* Lorentz factor *)
energycutoff=10; (* The maximum energy of the particle to which we limit the simulation *)
Print["Loading basic data"];

(* Monte Carlo simulation *) 
(* Necessary functions determining the particle's motion *)
Yscat[x_,{\[Epsilon]r_,x0_,e_,\[Lambda]_}]:=
	Module[{y,g2,g3,ys,y1,y2,y3,k2,X,X0},
		g2=1/12-1/\[Lambda]^2;
		g3=1/6^3-1/(12 \[Lambda]^2)-(e^2-1)/(4 \[Lambda]^2);
		ys=Sort[N[y/.Solve[4 y^3 - g2 y -g3==0,y]]];
		y1=ys[[1]];
		y2=ys[[2]];
		y3=ys[[3]];
		k2=(y2-y1)/(y3-y1);
		-\[Epsilon]r/Sqrt[y3-y1]*(EllipticF[ArcCos[Sqrt[(y2+1/12-1/(2x))/(y2-y1)]],k2]-EllipticF[ArcCos[Sqrt[(y2+1/12-1/(2x0))/(y2-y1)]],k2])
	]
	
Yabs[x_,{x0_,e_,\[Lambda]_}]:=
	Module[{y,g2,g3,poly,yreal,part2,q,p,mu,k2,X},
		g2=1/12-1/\[Lambda]^2;
		g3=1/6^3-1/(12 \[Lambda]^2)-(e^2-1)/(4 \[Lambda]^2);
		poly=4y^3-g2 y - g3;
		yreal=First[y/.Solve[poly==0,y,Reals,WorkingPrecision->10]];
		part2=FactorList[poly][[3,1]];
		q=SeriesCoefficient[part2,{y,0,0}];
		p=SeriesCoefficient[part2,{y,0,1}];
		mu=Sqrt[part2/.{y->yreal}];
		k2=(1/2)(1-(yreal+p/2)/mu);
		1/(2Sqrt[mu])*(EllipticF[2ArcTan[Sqrt[(-1/12+1/(2x)-yreal)/mu]],k2]-EllipticF[2ArcTan[Sqrt[(-1/12+1/(2x0)-yreal)/mu]],k2])
	]
(* The pericenter radius for unbound scattered orbit *)	
\[Xi]peri[\[Epsilon]_,\[Lambda]_]:=
	Module[{a0,a1,a2,a3,a4,f,\[Xi]},
		a0=(\[Epsilon]^2-1)/\[Lambda]^2;
		a1=1/(2 \[Lambda]^2);
		a2=-1/6;
		a3=1/2;
		a4=0;
		f=a0 \[Xi]^4 + 4 a1 \[Xi]^3 + 6 a2 \[Xi]^2 + 4 a3 \[Xi]+a4;
		Max[\[Xi]/.Solve[f==0,\[Xi]]]
	]
	
\[Lambda]c[\[Epsilon]_]:=Sqrt[12/(1-4/(1+3\[Epsilon]/Sqrt[9 \[Epsilon]^2 - 8])^2)] (* Critical angular momentum for a given energy \[Epsilon] *)
\[Lambda]max[\[Epsilon]_,\[Xi]_]:=\[Xi] Sqrt[\[Epsilon]^2 \[Xi]/(\[Xi]-2)-1] (* Maximal angular momentum for a geodesic of energy \[Epsilon] at the radius \[Xi] *)
\[Lambda]cutoff=\[Lambda]max[energycutoff,x\[Xi]0]; (* The maximum angular momentum of the particle to which we limit the simulation *)
Print["Loading the functions"];

(* Volume setup *) 

voltotM=NIntegrate[Exp[-\[Beta]/Sqrt[1-v^2]*(\[CurlyEpsilon]+v*Sqrt[\[CurlyEpsilon]^2-1]*Cos[\[Phi]])]\[Lambda]c[\[CurlyEpsilon]],{\[Phi],0,2*Pi},{\[CurlyEpsilon],1,energycutoff}];
voltotM2=NIntegrate[Exp[-\[Beta]/Sqrt[1-v^2]*(\[CurlyEpsilon]+v*Sqrt[\[CurlyEpsilon]^2-1]*Cos[\[Phi]])](\[Lambda]max[\[CurlyEpsilon],x\[Xi]0]-\[Lambda]c[\[CurlyEpsilon]]),{\[Phi],0,2*Pi},{\[CurlyEpsilon],1,energycutoff}];
voltotP=NIntegrate[Exp[-\[Beta]/Sqrt[1-v^2]*(\[CurlyEpsilon]-v*Sqrt[\[CurlyEpsilon]^2-1]*Cos[\[Phi]])]\[Lambda]c[\[CurlyEpsilon]],{\[Phi],0,2*Pi},{\[CurlyEpsilon],1,energycutoff}];
voltotP2=NIntegrate[Exp[-\[Beta]/Sqrt[1-v^2]*(\[CurlyEpsilon]-v*Sqrt[\[CurlyEpsilon]^2-1]*Cos[\[Phi]])](\[Lambda]max[\[CurlyEpsilon],x\[Xi]0]-\[Lambda]c[\[CurlyEpsilon]]),{\[Phi],0,2*Pi},{\[CurlyEpsilon],1,energycutoff}];
Print["Parallelization of definied functions"];

(* Couting trajectories of absorbed trajectories in a circle with constant radius r0 *) 

n\[Phi]=360; (* Number of cells around the BH *) 
nr=100; (* Number of districts in which the count was carried out *)
d\[Phi]=2Pi/n\[Phi]; (* Pixel width *) 
dr= (20-2)/nr; (* Pixel height *)
loopTime=SessionTime[]; (* Beginning of the time countdown to analyze the duration of individual loops *)
numberOfParts=20; (* Number of data files to read *)
nSubsets=$ProcessorCount; (* Number of subsets/processors - It can be adjusted manually or automatically by the $ProcessorCount command, then it takes all available ones*)

(* We import information about the total length of the data, without division into parts, to determine normalization *)
LengthScatP=Import[FileNameJoin[{notebookDir, "Data","dataLenghts.m"}]][[3]];
normalizScatP=LengthScatP*d\[Phi]; 

StepTime0=SessionTime[];

Print["r0=2.00001, time (in minutes):",(StepTime0-loopTime)/60];
licznik=1;

For[r0=2.00001, r0<=20.1, r0=r0+dr,
	JtScatResultsP=Table[0,{pixel,1,n\[Phi]}];
	JrScatResultsP=Table[0,{pixel,1,n\[Phi]}];
	JphiScatResultsP=Table[0,{pixel,1,n\[Phi]}];
	
	For[j=1,j<=numberOfParts,j++,
		(* We parallelize the function definitions again to make sure that they are loaded on each process *)
		LaunchKernels[];$KernelCount;
		DistributeDefinitions[Yscat, \[Xi]peri, n\[Phi], voltotM, voltotM2, voltotP, voltotP2, n\[Phi], numberOfParts, d\[Phi], dr];
		Print["Another parallelization of definitions"];
		
		(* We create the path to the data files *)
		fileNameScatP=FileNameJoin[{notebookDir,"Data", StringJoin["dataScatP_Part",ToString[j],".m"]}];
		(* Reading data from files *)
		scatP=Import[fileNameScatP,"Package"];
		(* Split 'data' into subsets. To ensure that all data is included, you can use the UpTo option in the Partition function to ensure that all elements are included in the subsets, even if the last subset has fewer elements than the rest.*)
		ScatPSubsets=Partition[scatP,UpTo[Ceiling[Length[scatP]/nSubsets]]];
		
		(* A function for processing a single subset of data *)
		processSubsetScatP[subset_]:=
		Module[{countsLocal=Table[0,{pixel,1,n\[Phi]}],countsrLocal=Table[0,{pixel,1,n\[Phi]}],countsphiLocal=Table[0,{pixel,1,n\[Phi]}]},
		
			For[i=1,i<=Length[subset],i++,
				(* Loading orbit parameters *)
				dir= subset[[i,4]];
				\[Lambda]=subset[[i,2]];
				\[Phi]0=subset[[i,3]];
				\[Epsilon]=subset[[i,1]];
				per=\[Xi]peri[\[Epsilon],\[Lambda]]//Re;
				(* We check whether the counting point is before the pericenter is located *)
				If[r0>per,
					(* We count the angle at which the particle meets the target *)
					psiIN=\[Phi]0-dir*Yscat[r0,{+1,x\[Xi]0,\[Epsilon],\[Lambda]}]//Re;
					psiIN2=Mod[psiIN,2Pi];
					pixelin=IntegerPart[psiIN2/d\[Phi]]+1;
					(* Counting the values of observables in a given pixel *)
					weight=\[Epsilon]/(r0*Sqrt[\[Epsilon]^2-(1-2/r0)(1+\[Lambda]^2/r0^2)]);
					weightr=1/(r0*(1-2/r0)); 
					weightphi=dir *\[Lambda]/(r0*Sqrt[\[Epsilon]^2-(1-2/r0)(1+\[Lambda]^2/r0^2)]);
					weight=2 weight* voltotP2/normalizScatP;
					weightrIN=-2  weightr*voltotP2/normalizScatP;
					weightphi=-2 weightphi*voltotP2/normalizScatP;
					countsLocal[[pixelin]]+=N[weight];
					countsrLocal[[pixelin]]+=N[weightrIN];
					countsphiLocal[[pixelin]]+=N[weightphi];
				];
			];
		{countsLocal,countsrLocal,countsphiLocal} (*Returning results as a list of lists*)
		];
		(*Parallel processing of each subset*)
		allResultsScatP=ParallelMap[processSubsetScatP,ScatPSubsets];
		(*Combining results into individual lists with data*)
		finalCountsScatP=N[Total[allResultsScatP[[All,1]]]];
		finalCountsrScatP=N[Total[allResultsScatP[[All,2]]]];
		finalCountsphiScatP=N[Total[allResultsScatP[[All,3]]]];
		
		(*Accumulation of results for each loop run*)
		JtScatResultsP=JtScatResultsP+finalCountsScatP;
		JrScatResultsP=JrScatResultsP+finalCountsrScatP;
		JphiScatResultsP=JphiScatResultsP+finalCountsphiScatP;
		
		(* Memory clearing *) 
		Clear[scatP, ScatPSubsets];
		ParallelEvaluate[Clear[scatP, ScatPSubsets, dir, \[Lambda], \[Phi]0, \[Epsilon], per, psiIN, psiIN2, pixelin, weight, weightr, weightphi, weight, weightrIN, weightphi, allResultsScatP, finalCountsScatP, finalCountsrScatP, finalCountsphiScatP ]];
		CloseKernels[];
		ClearSystemCache[];
		
		StepTime1=SessionTime[];
		Print["Number of finished inner loop:",j ", time (in minutes):",(StepTime1-StepTime0)/60];
		
	];

	(* Assigning the final values \:200b\:200bof currents J their angular position *)
	AllResultsJtScatP={};
	AllResultsJrScatP={};
	AllResultsJphiScatP={};
	
	For[pixel=1,pixel<=n\[Phi],pixel++,
		\[Phi]dn=(pixel-1)d\[Phi];
		\[Phi]up=pixel d\[Phi];
		\[Phi]sr=(\[Phi]dn+\[Phi]up)/2;
		AppendTo[AllResultsJtScatP,{N[\[Phi]sr], JtScatResultsP[[pixel]]}];
		AppendTo[AllResultsJrScatP,{N[\[Phi]sr], JrScatResultsP[[pixel]]}];
		AppendTo[AllResultsJphiScatP,{N[\[Phi]sr],JphiScatResultsP[[pixel]]}];
		];
	
	(* Exporting data to external files *)
	Export[FileNameJoin[{notebookDir, "Results", StringJoin["JtScatP",ToString[licznik],".m"]}],{r0,AllResultsJtScatP},"Package"];
	Export[FileNameJoin[{notebookDir, "Results", StringJoin["JrScatP",ToString[licznik],".m"]}],{r0,AllResultsJrScatP},"Package"];
	Export[FileNameJoin[{notebookDir, "Results", StringJoin["JphiScatP",ToString[licznik],".m"]}],{r0,AllResultsJphiScatP},"Package"];
	Print["Results file number:", licznik];
	
	licznik=licznik+1;
	
	StepTime2=SessionTime[];
	Print["r0=",r0+dr,", time (in minutes):",(StepTime2-StepTime0)/60];

];

voverallTime=(SessionTime[]-beginTime)/60
Print["the end and total time (in minutes):", voverallTime];

Quit[];
