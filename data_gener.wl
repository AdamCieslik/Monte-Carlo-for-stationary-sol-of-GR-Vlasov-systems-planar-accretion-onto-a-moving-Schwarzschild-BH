(* ::Package:: *)

(* Setup of the basic quantities *)
ClearAll["Global`*"]
beginTime=SessionTime[];
notebookDir="directory of the program file";
\[Xi]h=2; (* Location of the horizon *) 
\[Xi]f=3; (* Location of the photon orbit *) 
(* Mass of the test particle and the black hole *)
M = 1; 
m = 1;
(* \[Alpha] and \[Beta] from Maxwell-Juttner distribytion function *)
\[Alpha] = 1;
\[Beta] = 1; 

x\[Xi]0=1000; (* Outer radius of the numerical grid *)

v=0.95; (* Velocity of the boost *)
\[Gamma]=1/Sqrt[1-v^2]; (* Lorentz factor *)
energycutoff=10; (* The maximum energy of the particle to which we limit the simulation *)
Print["Loading basic data"];

(* Monte Carlo simulation *) 

(* Number of particles *)
num=2*10^8;

\[Lambda]c[\[Epsilon]_]:=Sqrt[12/(1-4/(1+3\[Epsilon]/Sqrt[9 \[Epsilon]^2 - 8])^2)] (* Critical angular momentum for a given energy \[Epsilon] *)
\[Lambda]max[\[Epsilon]_,\[Xi]_]:=\[Xi] Sqrt[\[Epsilon]^2 \[Xi]/(\[Xi]-2)-1] (* Maximal angular momentum for a geodesic of energy \[Epsilon] at the radius \[Xi] *)
\[Lambda]cutoff=\[Lambda]max[energycutoff,x\[Xi]0]; (* The maximum angular momentum of the particle to which we limit the simulation *)

Print["Loading the basic quantities"];

(* Data generation *) 
(* Calculus for \[Epsilon]r=-1 *)
normM=NIntegrate[Exp[-\[Beta]/Sqrt[1-v^2]*(\[CurlyEpsilon]+v*Sqrt[\[CurlyEpsilon]^2-1]*Cos[\[Phi]])],{\[Phi],0,2*Pi},{\[CurlyEpsilon],1,energycutoff}];
distM=ProbabilityDistribution[Exp[-\[Beta]/Sqrt[1-v^2]*(\[CurlyEpsilon]+v*Sqrt[\[CurlyEpsilon]^2-1]*Cos[\[Phi]])]/normM,{\[CurlyEpsilon],1,energycutoff},{\[Phi],0,2*Pi}];
(* Generating initial data *)
daneM0=.
(* Generating data from distribution *)
daneM0=RandomVariate[distM,num,Method->{"MCMC","Thinning"->2,"InitialVariance"->1}];
daneM=Transpose[{daneM0[[All,1]],RandomReal[{0,\[Lambda]cutoff},num],daneM0[[All,2]],RandomChoice[{-1,1},num]}];
Print["DistM data generation, time (in minutes):", (SessionTime[]-beginTime)/60]
(* Dividing data into scattered and absorbed parts *)
Clear[results,\[Epsilon],\[Lambda],\[Phi]0,dir]

results=Table[
	{\[Epsilon],\[Lambda],\[Phi]0,dir}=daneM[[i]][[1;;4]];
	\[Lambda]up=\[Lambda]max[\[Epsilon],x\[Xi]0];
	If[
		\[Lambda]<\[Lambda]max[\[Epsilon],x\[Xi]0],
		{\[Epsilon],\[Lambda],\[Phi]0,dir},Nothing
	],
	{i,1,num}];
(* Filtering the results *)
scatM=Select[results,#[[2]]>\[Lambda]c[#[[1]]]&];
absM=Select[results,#[[2]]<=\[Lambda]c[#[[1]]]&];
Print["Data division, time (minutes):", (SessionTime[]-beginTime)/60];
Print["Length of absM data", Length[absM]];
Print["Length of scatM data", Length[scatM]];

(* You need to divide the data into smaller files because if they are too large, Mathematica will have problems working with them *)

numberOfParts=20; (* Number of divisions *)
partLengthScatM=Ceiling[Length[scatM]/numberOfParts];(*Length of one part for scatM*)

(* Export data for ScatM to separate files *)
Do[
Export[FileNameJoin[{notebookDir,StringJoin["dataScatM_Part",ToString[i],".m"]}],scatM[[1+(i-1)*partLengthScatM;;Min[i*partLengthScatM,Length[scatM]]]],"Package"]
,{i,1,numberOfParts}
]

partLengthAbsM=Ceiling[Length[absM]/numberOfParts];(* Length of one part for absM *)

(* Export data for absM to separate files *)
Do[
Export[FileNameJoin[{notebookDir,StringJoin["dataAbsM_Part",ToString[i],".m"]}],absM[[1+(i-1)*partLengthAbsM;;Min[i*partLengthAbsM,Length[absM]]]],"Package"]
,{i,1,numberOfParts}
]
Print["AbsM and scatM data export, time (in minutes):", (SessionTime[]-beginTime)/60];
(* Calculus for \[Epsilon]r=-1 *)
normP=NIntegrate[Exp[-\[Beta]/Sqrt[1-v^2]*(\[CurlyEpsilon]-v*Sqrt[\[CurlyEpsilon]^2-1]*Cos[\[Phi]])],{\[Phi],0,2*Pi},{\[CurlyEpsilon],1,energycutoff}];
distP=ProbabilityDistribution[Exp[-\[Beta]/Sqrt[1-v^2]*(\[CurlyEpsilon]-v*Sqrt[\[CurlyEpsilon]^2-1]*Cos[\[Phi]])]/normP,{\[CurlyEpsilon],1,energycutoff},{\[Phi],0,2*Pi}];
(* Generating initial data *)
dataP0=.
(* Generating data from distribution *)
dataP0=RandomVariate[distP,num,Method->{"MCMC","Thinning"->2,"InitialVariance"->1}];
dataP=Transpose[{dataP0[[All,1]],RandomReal[{0,\[Lambda]cutoff},num],dataP0[[All,2]],RandomChoice[{-1,1},num]}];
Print["DistP data generation, time (in minutes):", (SessionTime[]-beginTime)/60];

(* Dividing data into scattered and absorbed parts *)
Clear[results,\[Epsilon],\[Lambda],\[Phi]0,dir]

results=Table[{\[Epsilon],\[Lambda],\[Phi]0,dir}=dataP[[i]][[1;;4]];
\[Lambda]up=\[Lambda]max[\[Epsilon],x\[Xi]0];
If[\[Lambda]<\[Lambda]max[\[Epsilon],x\[Xi]0],
{\[Epsilon],\[Lambda],\[Phi]0,dir},Nothing],{i,1,num}];
(* Filtering the results *)
scatP=Select[results,#[[2]]>\[Lambda]c[#[[1]]]&];
absP=Select[results,#[[2]]<=\[Lambda]c[#[[1]]]&];
Print["Amount of absP data", Length[absP]];
Print["Amount of scatP data", Length[scatP]];

(* Length of one part for scatP *)
partLengthScatP=Ceiling[Length[scatP]/numberOfParts];

(* Export data for scatP to separate files *)
Do[
Export[FileNameJoin[{notebookDir,StringJoin["dataScatP_Part",ToString[i],".m"]}],scatP[[1+(i-1)*partLengthScatP;;Min[i*partLengthScatP,Length[scatP]]]],"Package"]
,{i,1,numberOfParts}
]

(* Exporting information about the total length of files. It is needed to set normalization in the counting process. *)
Export[FileNameJoin[{notebookDir,"dataLenghts.m"}],{Length[absM],Length[scatM],Length[scatP]}]

Print["Exporting scatP data and information about the lengths of all data"];

voverallTime=(SessionTime[]-beginTime)/60
Print["The end and total time (in minutes):", voverallTime]

Quit[];
