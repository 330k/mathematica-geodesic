(* ::Package:: *)

(* :Mathematica Version: 8.0 *)

(* :Name: Geodesic`GeoGraphicLib` *)

(* :Title: GeographicLib Interface Package for Mathematica *)

(* :Author: Kei Misawa *)

(* :Summary:
 GeographicLib is written by Charles F.F. Karney.
 

 This software is released under the MIT License.
  http://opensource.org/licenses/mit-license.php
*)

(* :Context: Geodesic`GeoGraphicLib` *)

(* :Package Version: 1.0 *)

(* :Copyright: Copyright 2014, Kei Misawa
   Copyright 2013, Charles F.F. Karney
 *)

(* :Keywords: Geodesy *)

(* :Warning:  *)

(* :Sources:  *)


BeginPackage["Geodesic`GeographicLib`","JLink`"];


GeoDestinationGL::usage="";


GeoDistanceGL::usage="";


GeoDirectionGL::usage="";


Begin["`Private`"];


InstallJava[];
AddToClassPath[FileNameJoin[{DirectoryName[FindFile["Geodesic`GeographicLib`"]],"Java"}]];
LoadJavaClass["net.sf.geographiclib.Geodesic",AllowShortContext->False];
$wgs84=net`sf`geographiclib`Geodesic`WGS84;


geoDirect[{lat_,lon_},{d_,\[Lambda]_},datum_,opts:OptionsPattern[]]:=JavaBlock[Module[{
	a,b,f,geodesic,result
},
	{a,b,f}=ellipsoidParameters[datum];
	geodesic=JavaNew["net.sf.geographiclib.Geodesic",a,f];
	result=geodesic@Direct[lat,lon,\[Lambda],d];
	GeoPosition[{result@lat2,result@lon2},datum]
]];


geoInverse[{lat1_,lon1_},{lat2_,lon2_},datum_,opts:OptionsPattern[]]:=JavaBlock[Module[{
	a,b,f,geodesic,result
},
	{a,b,f}=ellipsoidParameters[datum];
	geodesic=JavaNew["net.sf.geographiclib.Geodesic",a,f];
	result=geodesic@Inverse[lat1,lon1,lat2,lon2];
	{result@s12,result@azi1}
]];


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
Options[GeoDestinationGL]=Options[geoDirect];
GeoDestinationGL[{lat_?validLatQ,lon_?validLonQ,h___},{d_?NumericQ/;0<=d,\[Lambda]_?NumericQ},opts:OptionsPattern[]]:=geoDirect[{lat,lon},{d,\[Lambda]},"ITRF00",opts];
GeoDestinationGL[HoldPattern[GeoPosition[{lat_?validLatQ,lon_?validLonQ,h___},datum_:"ITRF00"]],{d_?NumericQ/;0<=d,\[Lambda]_?NumericQ},opts:OptionsPattern[]]:=geoDirect[{lat,lon},{d,\[Lambda]},datum,opts];
SyntaxInformation[GeoDestinationGL]={"ArgumentsPattern"->{_,{_,_},OptionsPattern[]}};


Options[GeoDistanceGL]=Options[geoInverse];
GeoDistanceGL[{lat1_?validLatQ,lon1_?validLonQ,h___},{lat2_?validLatQ,lon2_?validLonQ,h___},datum_?datumQ,opts:OptionsPattern[]]:=geoInverse[{lat1,lon1},{lat2,lon2},datum,opts][[1]];
GeoDistanceGL[{lat1_?validLatQ,lon1_?validLonQ,h___},{lat2_?validLatQ,lon2_?validLonQ,h___},opts:OptionsPattern[]]:=geoInverse[{lat1,lon1},{lat2,lon2},"ITRF00",opts][[1]];
GeoDistanceGL[
	HoldPattern[GeoPosition[{lat1_?validLatQ,lon1_?validLonQ,h1___},datum_:"ITRF00"]],
	HoldPattern[GeoPosition[{lat2_?validLatQ,lon2_?validLonQ,h2___},datum_:"ITRF00"]],
	opts:OptionsPattern[]
]:=geoInverse[{lat1,lon1},{lat2,lon2},datum,opts][[1]];
SyntaxInformation[GeoDistanceGL]={"ArgumentsPattern"->{{_,_},{_,_},_.,OptionsPattern[]}};


Options[GeoDirectionGL]=Options[geoInverse];
GeoDirectionGL[{lat1_?validLatQ,lon1_?validLonQ,h___},{lat2_?validLatQ,lon2_?validLonQ,h___},datum_?datumQ,opts:OptionsPattern[]]:=geoInverse[{lat1,lon1},{lat2,lon2},datum,opts][[2]];
GeoDirectionGL[{lat1_?validLatQ,lon1_?validLonQ,h___},{lat2_?validLatQ,lon2_?validLonQ,h___},opts:OptionsPattern[]]:=geoInverse[{lat1,lon1},{lat2,lon2},"ITRF00",opts][[2]];
GeoDirectionGL[
	HoldPattern[GeoPosition[{lat1_?validLatQ,lon1_?validLonQ,h1___},datum_:"ITRF00"]],
	HoldPattern[GeoPosition[{lat2_?validLatQ,lon2_?validLonQ,h2___},datum_:"ITRF00"]],
	opts:OptionsPattern[]
]:=geoInverse[{lat1,lon1},{lat2,lon2},datum,opts][[2]];
SyntaxInformation[GeoDirectionGL]={"ArgumentsPattern"->{{_,_},{_,_},_.,OptionsPattern[]}};


End[];


EndPackage[];
