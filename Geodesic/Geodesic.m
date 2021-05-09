(* ::Package:: *)

(* :Mathematica Version: 8.0 *)

(* :Name: Geodesic` *)

(* :Title: Alternative Geodesic Calculation Package *)

(* :Author: Kei Misawa *)

(* :Summary:
  This package provides 3 functions, GeoDestination2, GeoDistance2 and
 GeoDirection2, which are substitute for GeoDestination, GeoDistance and
 GeoDirection respectively. This package can solve geodesic problems
 more accurately with arbitrary precision calculation by
 designating WorkingPrecision option.

  Although this package and GeographicLib refer to the same paper by
 Charles F.F.Karney[1,2], they differ in the way to calculate.
 GeographicLib uses functional expansion to calculate elliptic
 integrals, but this package uses Mathematica built-in functions such as 
 EllipticE and EllipticF. Additionally, inverse function of elliptic
 integral is represented by InverseFunction, and FindRoot is used to
 solve the inverse problems. So this "make Mathematica calculate" concept
 enables arbitrary precision calculation without numerical techniques, but
 the precision of calculation depends on Mathematica.

  GeoDestination2, GeoDistance2 and GeoDirection2 can be used almost
 samely as GeoDestination, GeoDistance and GeoDirection.
 However, note that not all patterns of arguments are supported.
 For example, GeoPositionXYZ and GeoPositionENU are not supported.

  In order to calculate with arbitrary precision, set option value to
 WorkingPrecision.

 This software is released under the MIT License.
  http://opensource.org/licenses/mit-license.php
*)

(* :Context: Geodesic` *)

(* :Package Version: 1.0 *)

(* :Copyright: Copyright 2014, Kei Misawa *)

(* :Keywords: Geodesy *)

(* :Warning:  *)

(* :Sources: 
  [1] Charles F.F. Karney, "Algorithms for geodesics," J. Geodesy 87, 43-55 (2013),
      http://dx.doi.org/10.1007/s00190-012-0578-z
  [2] Charles F.F. Karney, "Geodesics on an ellipsoid of revolution," Feb. 2011,
      http://arxiv.org/abs/1102.1215v1
*)


BeginPackage["Geodesic`"];


GeoDestination2::usage="GeoDestination2[\!\(\*
StyleBox[\"pos\", \"TI\"]\),{\!\(\*
StyleBox[\"d\", \"TI\"]\),\!\(\*
StyleBox[\"\[Alpha]\", \"TR\"]\)}] returns position from \!\(\*
StyleBox[\"pos\", \"TI\"]\) to direction \!\(\*
StyleBox[\"\[Alpha]\", \"TR\"]\) and distance \!\(\*
StyleBox[\"d\", \"TI\"]\).";


GeoDistance2::usage="GeoDistance2[{\!\(\*SubscriptBox[
StyleBox[\"lat\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\),\!\(\*SubscriptBox[
StyleBox[\"long\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\)},{\!\(\*SubscriptBox[
StyleBox[\"lat\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\),\!\(\*SubscriptBox[
StyleBox[\"long\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\)}] returns geodetic distance between 2 points.\nOptions are equal to FindRoot.";


GeoDirection2::usage="GeoDirection2[{\!\(\*SubscriptBox[
StyleBox[\"lat\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\),\!\(\*SubscriptBox[
StyleBox[\"long\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\)},{\!\(\*SubscriptBox[
StyleBox[\"lat\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\),\!\(\*SubscriptBox[
StyleBox[\"long\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\)}] returns direction from 1st point to 2nd point.\nOptions are equal to FindRoot.";


Begin["`Private`"];


(* internal function which solves direct problems *)
Options[geoDirect]={
	WorkingPrecision->MachinePrecision
};
geoDirect[{lat_,lon_},{d_,\[Lambda]_},datum_,opts:OptionsPattern[]]:=Module[{
	a,f,b,e2,eprime2,\[Phi]1,\[Alpha]1,s12,\[Beta]1,
	\[Alpha]0,\[Sigma]1,\[Omega]1,k2,s1,s2,\[Sigma]2,\[Alpha]2,
	\[Beta]2,\[Omega]2,\[Lambda]1,\[Lambda]2,\[Lambda]12,\[Phi]2
},
	{a,b,f}=SetPrecision[ellipsoidParameters[datum],OptionValue[WorkingPrecision]];
	e2=f(2-f);
	eprime2=e2/(1-e2);

	\[Phi]1=SetPrecision[lat Degree,OptionValue[WorkingPrecision]];
	\[Alpha]1=SetPrecision[Mod[\[Lambda],360]Degree,OptionValue[WorkingPrecision]];
	\[Lambda]1=SetPrecision[lon Degree,OptionValue[WorkingPrecision]];
	s12=SetPrecision[d,OptionValue[WorkingPrecision]];

	\[Beta]1=ArcTan[(1-f)Tan[\[Phi]1]];
	\[Alpha]0=ArcSin[Sin[\[Alpha]1]Cos[\[Beta]1]];
	\[Sigma]1=ArcTan[Cos[\[Alpha]1]Cos[\[Beta]1],Sin[\[Beta]1]];
	\[Omega]1=ArcTan[Cos[\[Sigma]1],Sin[\[Alpha]0]Sin[\[Sigma]1]];
	k2=eprime2 Cos[\[Alpha]0]^2;
	s1=b Re@EllipticE[\[Sigma]1,-k2];
	s2=s1+s12;
	\[Sigma]2=InverseFunction[EllipticE,1,2][s2/b,-k2];
	\[Alpha]2=ArcTan[Cos[\[Alpha]0]Cos[\[Sigma]2],Sin[\[Alpha]0]];
	\[Beta]2=ArcTan[Sqrt[(Cos[\[Alpha]0]Cos[\[Sigma]2])^2+Sin[\[Alpha]0]^2],Cos[\[Alpha]0]Sin[\[Sigma]2]];
	\[Omega]2=ArcTan[Cos[\[Sigma]2],Sin[\[Alpha]0]Sin[\[Sigma]2]];
	\[Lambda]12=(1-f)Sin[\[Alpha]0](G[\[Sigma]2,Cos[\[Alpha]0]^2,-k2]-G[\[Sigma]1,Cos[\[Alpha]0]^2,-k2]);
	\[Phi]2=ArcTan[Tan[\[Beta]2]/(1-f)];
	\[Lambda]2=Mod[\[Lambda]1+\[Lambda]12,2Pi];

	GeoPosition[{\[Phi]2/Degree,\[Lambda]2/Degree},datum]
];


(* internal function which solves inverse problems *)
Options[geoInverse]=Options[FindRoot];(*{
	WorkingPrecision->MachinePrecision,
	AccuracyGoal->Automatic,
	PrecisionGoal->Automatic,
	MaxIterations->100
};*)
geoInverse[{lat1_,lon1_},{lat2_,lon2_},datum_,opts:OptionsPattern[]]:=Catch[Module[{
	a,b,f,e2,eprime2,inverselat=False,inverselon=False,swapped=False,iterations=0,
	\[Phi]1,\[Phi]2,\[Lambda]12org,\[Beta]1,\[Beta]2,\[Alpha]1,\[Omega]12,\[Alpha]10,\[Lambda]12,
	\[Alpha]0,\[Sigma]1,\[Omega]1,k2,\[Alpha]2,\[Sigma]2,\[Omega]2,
	s1,s2,s12
},
	{a,b,f}=SetPrecision[ellipsoidParameters[datum],OptionValue[WorkingPrecision]];
	e2=SetPrecision[f(2-f),OptionValue[WorkingPrecision]];
	eprime2=SetPrecision[e2/(1-e2),OptionValue[WorkingPrecision]];

	\[Phi]1=SetPrecision[lat1 Degree,OptionValue[WorkingPrecision]];
	\[Phi]2=SetPrecision[lat2 Degree,OptionValue[WorkingPrecision]];
	\[Lambda]12org=SetPrecision[(Mod[lon2+180,360]-Mod[lon1+180,360])Degree,OptionValue[WorkingPrecision]];
	If[TrueQ[Abs[\[Phi]1]<Abs[\[Phi]2]],
		{\[Phi]1,\[Phi]2}={\[Phi]2,\[Phi]1};
		swapped=True;
	];
	If[TrueQ[\[Phi]1>0],
		\[Phi]1=-\[Phi]1;
		\[Phi]2=-\[Phi]2;
		inverselat=True;
	];
	If[TrueQ[lon2-lon1<0],
		\[Lambda]12org=-\[Lambda]12org;
		inverselon=True;
	];
	If[TrueQ[\[Lambda]12org>Pi],
		\[Lambda]12org=2Pi-\[Lambda]12org;
		inverselon=!inverselon;
	];
	\[Beta]1=ArcTan[(1-f)Tan[\[Phi]1]];
	\[Beta]2=ArcTan[(1-f)Tan[\[Phi]2]];

	{s12,\[Alpha]1,\[Alpha]2}=Catch[
		Which[
			TrueQ[\[Beta]1==\[Beta]2==0&&Abs[\[Lambda]12org]<(1-f)*180*Degree],
			(* equator *)
			(*\[Alpha]1=Pi/2;
			\[Alpha]2=Pi/2;
			s12=a Abs[\[Lambda]12org];Print[FullForm@s12];*)
			Throw[{a*Abs[\[Lambda]12org],Pi/2,Pi/2}];
			,
			TrueQ[\[Beta]1==-Pi/2||\[Lambda]12org==-Pi],
			(* meridian *)
			\[Alpha]1=Pi;
			\[Alpha]2=0;
			\[Sigma]1=ArcTan[Cos[\[Alpha]1]Cos[\[Beta]1],Sin[\[Beta]1]];
			\[Sigma]2=ArcTan[Cos[\[Alpha]2]Cos[\[Beta]2],Sin[\[Beta]2]];
			k2=eprime2;
			s1=b Re[EllipticE[\[Sigma]1,-k2]];
			s2=b Re[EllipticE[\[Sigma]2,-k2]];
			s12=s2-s1;
			Throw[{s2-s1,Pi,0}];
		];
		Which[
			TrueQ[Abs[\[Beta]1-\[Beta]2(*Abs[\[Beta]1]-Abs[\[Beta]2]*)]<f*Degree&&Abs[\[Lambda]12org]>(1-f)*180*Degree],
			(* nearly antipodal *)
			Module[{\[CapitalDelta],x,y,\[Mu]},
				\[CapitalDelta]=f a Pi Cos[\[Beta]1]^2;
				x=(\[Lambda]12org-Pi)a Cos[\[Beta]1]/\[CapitalDelta];
				y=(\[Beta]1+\[Beta]2)a/\[CapitalDelta];
				If[TrueQ[y!=0],
					Module[{\[Mu]ans},
						\[Mu]ans=NSolve[\[Mu]^4+2\[Mu]^3+(1-x^2-y^2)\[Mu]^2-2y^2\[Mu]-y^2==0&&0<=\[Mu],\[Mu]];
						If[Length[\[Mu]ans]>=1,
							\[Mu]=First[\[Mu]/.\[Mu]ans];
							\[Alpha]10=ArcTan[y/\[Mu],-x/(1+\[Mu])];
							,
							\[Alpha]10=ArcTan[Sqrt[Max[0,1-x^2]],-x];
						];
					];
					(*Print[{Abs[\[Beta]1-\[Beta]2],x,y,\[Mu],\[Alpha]10}];*)
					,
					\[Alpha]10=ArcTan[Sqrt[Max[0,1-x^2]],-x];
				]
			];
			,
			True,
			(* else *)
			\[Omega]12=\[Lambda]12org/Sqrt[1-e2 ((Cos[\[Beta]1]+Cos[\[Beta]2])/2)^2];
			\[Alpha]10=ArcTan[Cos[\[Beta]1]Sin[\[Beta]2]-Sin[\[Beta]1]Cos[\[Beta]2]Cos[\[Omega]12],Cos[\[Beta]2]Sin[\[Omega]12]];
			(*Print[{\[Omega]12,\[Alpha]10,ArcTan[Cos[\[Beta]1]Sin[\[Beta]2]-Sin[\[Beta]1]Cos[\[Beta]2]Cos[\[Omega]12],Cos[\[Beta]2]Sin[\[Omega]12]]}];*)
		];

		(* Solve \[Alpha]1 *)
		(*\[Lambda]12=(-1+f) Cos[\[Beta]1] (G[ArcTan[Cos[\[Alpha]1] Cos[\[Beta]1],Sin[\[Beta]1]],1-Cos[\[Beta]1]^2 Sin[\[Alpha]1]^2,eprime2 (-1+Cos[\[Beta]1]^2 Sin[\[Alpha]1]^2)]-G[ArcTan[Sqrt[Cos[\[Beta]2]^2-Cos[\[Beta]1]^2 Sin[\[Alpha]1]^2],Sin[\[Beta]2]],1-Cos[\[Beta]1]^2 Sin[\[Alpha]1]^2,eprime2 (-1+Cos[\[Beta]1]^2 Sin[\[Alpha]1]^2)]) Sin[\[Alpha]1];*)
		\[Lambda]12=-(-1+f) Cos[\[Beta]1] (eprime2 EllipticF[ArcTan[Cos[\[Alpha]1] Cos[\[Beta]1],Sin[\[Beta]1]],eprime2 (-1+Cos[\[Beta]1]^2 Sin[\[Alpha]1]^2)]-eprime2 EllipticF[ArcTan[Sqrt[Cos[\[Beta]2]^2-Cos[\[Beta]1]^2 Sin[\[Alpha]1]^2],Sin[\[Beta]2]],eprime2 (-1+Cos[\[Beta]1]^2 Sin[\[Alpha]1]^2)]-(1+eprime2) (EllipticPi[1-Cos[\[Beta]1]^2 Sin[\[Alpha]1]^2,ArcTan[Cos[\[Alpha]1] Cos[\[Beta]1],Sin[\[Beta]1]],eprime2 (-1+Cos[\[Beta]1]^2 Sin[\[Alpha]1]^2)]-EllipticPi[1-Cos[\[Beta]1]^2 Sin[\[Alpha]1]^2,ArcTan[Sqrt[Cos[\[Beta]2]^2-Cos[\[Beta]1]^2 Sin[\[Alpha]1]^2],Sin[\[Beta]2]],eprime2 (-1+Cos[\[Beta]1]^2 Sin[\[Alpha]1]^2)])) Sin[\[Alpha]1];
		If[TrueQ[\[Beta]1==0],
			\[Lambda]12=-\[Lambda]12;
			\[Alpha]10=Pi/2+0.001;
		];
		(*Print[{f,eprime2,\[Beta]1,\[Beta]2,\[Lambda]12}];*)
		(*Print[Plot[Evaluate@{\[Lambda]12,\[Lambda]12org},{\[Alpha]1,-Pi,Pi}(*,PlotRange->{\[Lambda]12org*0.95,\[Lambda]12org*1.05}*)]];(* dubug *)*)
		(*Print[\[Alpha]10];*)
		\[Lambda]12=SetPrecision[\[Lambda]12,OptionValue[WorkingPrecision]];
		\[Alpha]1=\[Alpha]1/.First[FindRoot[
			\[Lambda]12==\[Lambda]12org
			,{\[Alpha]1,\[Alpha]10},
			Evaluate@FilterRules[{
				(*StepMonitor:>Print[NumberForm[\[Alpha]1,16]],*)
				opts},Options[FindRoot]
			]
		]];
		If[Not@NumberQ[\[Alpha]1],Throw[{$Failed,$Failed,$Failed}]];
		\[Alpha]0=ArcSin[Sin[\[Alpha]1]Cos[\[Beta]1]];
		\[Sigma]1=ArcTan[Cos[\[Alpha]1]Cos[\[Beta]1],Sin[\[Beta]1]];
		If[TrueQ[\[Beta]1==0],
			(* modify when \[Beta]1 == 0 (equator and antipodal ) *)
			\[Sigma]1=-\[Sigma]1;
			\[Alpha]1=Pi-\[Alpha]1;
		];
		\[Sigma]2=ArcTan[Sqrt[Cos[\[Alpha]1]^2Cos[\[Beta]1]^2+(Cos[\[Beta]2]^2-Cos[\[Beta]1]^2)],Sin[\[Beta]2]];
		\[Alpha]2=ArcTan[Cos[\[Alpha]0]Cos[\[Sigma]2],Sin[\[Alpha]0]];
		k2=eprime2 Cos[\[Alpha]0]^2;
		s1=b Re[EllipticE[\[Sigma]1,-k2]];
		s2=b Re[EllipticE[\[Sigma]2,-k2]];
		s12=s2-s1;
		Throw[{s2-s1,\[Alpha]1,\[Alpha]2}];
	];

	If[swapped,{\[Alpha]1,\[Alpha]2}={Pi-\[Alpha]2,Pi-\[Alpha]1};];
	If[inverselat,{\[Alpha]1,\[Alpha]2}={Pi-\[Alpha]1,Pi-\[Alpha]2};];
	If[inverselon,{\[Alpha]1,\[Alpha]2}={-\[Alpha]1,-\[Alpha]2}];
	\[Alpha]1=Mod[\[Alpha]1,2Pi];
	\[Alpha]2=Mod[\[Alpha]2,2Pi];
	Throw[{s12,\[Alpha]1/Degree,\[Alpha]2/Degree}];
]];


(* internal elliptic integral function *)
G[\[Phi]_,n_,m_]:=m/n EllipticF[\[Phi],m]+(1-m/n)EllipticPi[n,\[Phi],m];
G[\[Phi]_,0,m_]:=Evaluate[Limit[G[\[Phi],n,m],n->0]];


(* internal functions which check arguments *)
datumnames=GeodesyData[];
datumQ[datumstring_String]:=MemberQ[datumnames,datumstring];
datumQ[{a_?NumericQ,b_?NumericQ}]:=True/;(a>0&&b>0)
datumQ[___]:=False;


validLatQ[lat_?NumericQ]:=TrueQ[Element[lat,Reals]]&&-90<=lat<=90;
validLatQ[___]:=False;


validLonQ[lon_?NumericQ]:=TrueQ[Element[lon,Reals]];
validLonQ[___]:=False;


(* internal function which get ellipsoid parameters from datum *)
(ellipsoidParameters[#]:={
	GeodesyData[#,"SemimajorAxis"],
	GeodesyData[#,"SemiminorAxis"],
	1/GeodesyData[#,"InverseFlattening"]
})&/@datumnames;
ellipsoidParameters[{a_?Positive,b_?Positive}]:={a,b,1-b/a};


(* exported functions *)
Options[GeoDestination2]=Options[geoDirect];
GeoDestination2[{lat_?validLatQ,lon_?validLonQ,h___},{d_?NumericQ/;0<=d,\[Lambda]_?NumericQ},opts:OptionsPattern[]]:=geoDirect[{lat,lon},{d,\[Lambda]},"ITRF00",opts];
GeoDestination2[HoldPattern[GeoPosition[{lat_?validLatQ,lon_?validLonQ,h___},datum_:"ITRF00"]],{d_?NumericQ/;0<=d,\[Lambda]_?NumericQ},opts:OptionsPattern[]]:=geoDirect[{lat,lon},{d,\[Lambda]},datum,opts];
SyntaxInformation[GeoDestination2]={"ArgumentsPattern"->{_,{_,_},OptionsPattern[]}};


Options[GeoDistance2]=Options[geoInverse];
GeoDistance2[{lat1_?validLatQ,lon1_?validLonQ,h___},{lat2_?validLatQ,lon2_?validLonQ,h___},datum_?datumQ,opts:OptionsPattern[]]:=geoInverse[{lat1,lon1},{lat2,lon2},datum,opts][[1]];
GeoDistance2[{lat1_?validLatQ,lon1_?validLonQ,h___},{lat2_?validLatQ,lon2_?validLonQ,h___},opts:OptionsPattern[]]:=geoInverse[{lat1,lon1},{lat2,lon2},"ITRF00",opts][[1]];
GeoDistance2[
	HoldPattern[GeoPosition[{lat1_?validLatQ,lon1_?validLonQ,h1___},datum_:"ITRF00"]],
	HoldPattern[GeoPosition[{lat2_?validLatQ,lon2_?validLonQ,h2___},datum_:"ITRF00"]],
	opts:OptionsPattern[]
]:=geoInverse[{lat1,lon1},{lat2,lon2},datum,opts][[1]];
SyntaxInformation[GeoDistance2]={"ArgumentsPattern"->{{_,_},{_,_},_.,OptionsPattern[]}};


Options[GeoDirection2]=Options[geoInverse];
GeoDirection2[{lat1_?validLatQ,lon1_?validLonQ,h___},{lat2_?validLatQ,lon2_?validLonQ,h___},datum_?datumQ,opts:OptionsPattern[]]:=geoInverse[{lat1,lon1},{lat2,lon2},datum,opts][[2]];
GeoDirection2[{lat1_?validLatQ,lon1_?validLonQ,h___},{lat2_?validLatQ,lon2_?validLonQ,h___},opts:OptionsPattern[]]:=geoInverse[{lat1,lon1},{lat2,lon2},"ITRF00",opts][[2]];
GeoDirection2[
	HoldPattern[GeoPosition[{lat1_?validLatQ,lon1_?validLonQ,h1___},datum_:"ITRF00"]],
	HoldPattern[GeoPosition[{lat2_?validLatQ,lon2_?validLonQ,h2___},datum_:"ITRF00"]],
	opts:OptionsPattern[]
]:=geoInverse[{lat1,lon1},{lat2,lon2},datum,opts][[2]];
SyntaxInformation[GeoDirection2]={"ArgumentsPattern"->{{_,_},{_,_},_.,OptionsPattern[]}};


End[];


EndPackage[];
