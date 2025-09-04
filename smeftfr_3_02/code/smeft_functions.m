(* SmeftFR v3.0 package *)
(* Auxiliary functions to simplify and transform SMEFT Lagrangian *)


NumberToString = Function[{x},
(* string display format for floating point numbers *)
Module[{tmp,tmp1,tmp2},

If[ NumberQ[x],
  If[ Head[x] === Real || Head[x] === Rational,
    tmp1 = MantissaExponent[N[x]];
    tmp  = ToString[tmp1[[1]]] <> "*^" <> ToString[tmp1[[2]]];
  ];
  If[ Head[x] === Integer,
    tmp1 = MantissaExponent[x];
    tmp  = ToString[tmp1[[1]]] <> "*^" <> ToString[tmp1[[2]]];
  ];
  If[ Head[x] === Complex,
    tmp1 = MantissaExponent[N[Re[x]]];
    tmp1 = ToString[tmp1[[1]]] <> "*^" <> ToString[tmp1[[2]]];
    tmp = If[ Im[x] < 0, " - ", " + "];
    tmp2 = MantissaExponent[N[Abs[Im[x]]]];
    tmp2 = ToString[tmp2[[1]]] <> "*^" <> ToString[tmp2[[2]]];
    tmp  = tmp1 <> tmp <> tmp2 <> "*I";
  ];
,
  Print[ "NumberToString function called with argument ",x];
  Print[ "Argument is not a number, please correct!"];
  Abort[];
];

tmp

]
(* end of NumberToString *)
]




SMEFTExpReduce = Function[{x},
(* cuts higher order terms in 1/Lambda *)
Module[{k}, x /. Lam^k_ -> If[ k>SMEFT$ExpansionOrder, 0, Lam^k ] ]
]



SMEFTExpOrder = Function[{x},
(* expands expression to given order in 1/Lambda *)
Module[{tmp,k},

tmp = x /. Lam^k_ -> If[ k>SMEFT$ExpansionOrder, 0, Lam^k ];
Normal[ Series[ tmp, {Lam, 0, SMEFT$ExpansionOrder}] ]

]
(* end of SMEFTExpOrder *)
]



MBScalarList = Function[{op},
(* returns scalar operator in mass basis if present in SMEFT$OperatorList or 0 otherwise *)
If[MemberQ[ SMEFT$OperatorList, ToString[op] ], ToExpression[SMEFT$MB<>ToString[op]], 0]
(*end of MBScalarList *)			
]


WBScalarList = Function[{op},
(* returns scalar operator in Warsaw basis if present in SMEFT$OperatorList or 0 otherwise *)
If[MemberQ[ SMEFT$OperatorList, ToString[op] ], ToExpression[SMEFT$WB<>ToString[op]], 0]
(*end of WBScalarList *)			
]

MBTensor2List = Function[{op,i,j},
(* returns flavor operator with 2 indices in mass basis if present in SMEFT$OperatorList or 0 otherwise *)
If[MemberQ[ SMEFT$OperatorList, ToString[op] ],
            ToExpression[SMEFT$MB<>ToString[op]<>"["<>ToString[i]<>","<>ToString[j]<>"]"], 0]
(*end of MBTensor2List *)			
]

WBTensor2List = Function[{op,i,j},
(* returns flavor operator with 2 indices in Warsaw basis if present in SMEFT$OperatorList or 0 otherwise *)
If[MemberQ[ SMEFT$OperatorList, ToString[op] ],
            ToExpression[SMEFT$WB<>ToString[op]<>"["<>ToString[i]<>","<>ToString[j]<>"]"], 0]
(*end of WBTensor2List *)			
]


MBTensor4List = Function[{op,i,j,k,l},
(* returns flavor operator with 4 indices in mass basis if present in SMEFT$OperatorList or 0 otherwise *)
If[MemberQ[ SMEFT$OperatorList, ToString[op] ],
            ToExpression[SMEFT$MB<>ToString[op]<>"["<>ToString[i]<>","<>ToString[j]<>","<>ToString[k]<>","<>ToString[l]<>"]"], 0]
(*end of MBTensor4List *)			
]


WBTensor4List = Function[{op,i,j,k,l},
(* returns flavor operator with 4 indices in Warsaw basis if present in SMEFT$OperatorList
   or 0 otherwise *)

If[MemberQ[ SMEFT$OperatorList, ToString[op] ],
            ToExpression[SMEFT$WB<>ToString[op]<>"["<>ToString[i]<>","<>ToString[j]<>","<>ToString[k]<>","<>ToString[l]<>"]"], 0]

(*end of WBTensor4List *)			
]



GenerateOperatorLists[ OptionsPattern[ { Letter -> SMEFT$MB } ] ] :=
(* sublists of active operators *)
Module[{l},

l = OptionValue[Letter];
       
SMEFT$Dim6BosonOperators = ToExpression[l <> #] & /@ Intersection[
    Join[ GaugeBilinearOperators, GaugeTripleOperators ], SMEFT$OperatorList ];

SMEFT$Dim6FermionOperators = ToExpression[l <> #] & /@ Intersection[
    Join[ TwoFermionMassOperators, TwoFermionOperators,
          FourFermionOperators, BLViolatingOperators ], SMEFT$OperatorList ];

SMEFT$Dim8BosonOperators = ToExpression[l <> #] & /@ Intersection[
    Join[ GaugeBilinearOperators8, GaugeQuadrupleOperators ], SMEFT$OperatorList ];

SMEFT$Dim6NullList = Join[ (# -> 0 & /@ SMEFT$Dim6BosonOperators),
                           (#[__] -> 0 & /@ SMEFT$Dim6FermionOperators) ];

SMEFT$Dim8NullList = (# -> 0 & /@ SMEFT$Dim8BosonOperators);

]



CombineCommonPowers = Function[{expr,par,n},
(* rearranges power expansion into a_0 + a_1 par^1 + ... a_n par^n *)
Module[{tmp,fact,aa},
tmp =  expr /. par -> 0 // Simplify; 
If[ n > 0,
  For[aa=1, aa < n+1, aa++,
    fact = Coefficient[expr, par^aa]; 
    tmp = tmp + par^aa If[ aa==1, Dim6CoefficientsSimplify[fact], Simplify[fact, TimeConstraint->5] ]; 
  ];
];

tmp

]
(* end of CombineCommonPowers *)
]



Dim6CoefficientsSimplify = Function[{expr},
(* simplifies coefficients of dim 6 operators *)
Module[{tmp,fact,aa},

tmp = 0;
For[aa=1, aa < Length[SMEFT$Dim6NullList]+1, aa++,
    tmp = tmp + Simplify[expr /. Drop[SMEFT$Dim6NullList,{aa} ], TimeConstraint->10 ];
  ];

tmp

]
(* end of Dim6CoefficientsSimplify *)
]



CouplingsFormSimplify = Function[{vlist},
(* Simplify couplings in front of interaction vertices *)
Module[{tmp},				

tmp = Map[ CombineCommonPowers[#, Lam, SMEFT$ExpansionOrder] &, vlist[[All, 2]]];
Partition[Riffle[vlist[[All,1]], tmp],2]

]
(* end of CouplingsFormSimplify *)
]				



SU2Simplify = Function[{expr},
(*** SU(2) generators T^a_bc T^a_de reduction ***)
Module[{tmp,aa,y,z,y1,z1},

expr //. T[aa_,y_,z_] T[aa_,y1_,z1_] ->
      1/2 ( IndexDelta[y,z1] IndexDelta[z,y1] -
      1/3 IndexDelta[y,z] IndexDelta[y1,z1] ) // Expand

]
(* end of SU2Simplify *)
];




SMEFTExpandVertices[ OptionsPattern[{ Input -> "smeft",
				      ExpOrder-> SMEFT$ExpansionOrder,
				      ExpName -> "Exp",
				      Silent -> True}] ] :=
Module[{input, order, CPUTime, a, b},
(* expand normalization constant in terms of chosen set of input parameters *)

input = OptionValue[Input];
If [ ! MemberQ[{"smeft","user","none"}, input],
     Print["Currently only \"none\", \"smeft\" or  \"user\" allowed as options of SMEFTExpandVertices"];
     Abort[];
   ];

If [ ! MemberQ[ {0,1,2}, OptionValue[ExpOrder] ],
     Print["Currently only option ExpOrder = 0 (SM), ExpOrder = 1 (dim6) or = 2 (dim8) in SMEFTExpandVertices allowed"];
     Abort[];
   ];

(* store general expansion order *)       
order = SMEFT$ExpansionOrder;

(* initialize time counter *)
CPUTime = TimeUsed[];

NormalizationConstantsExpansion[input];

If [ OptionValue[ExpName] == "Latex",
  SMEFT$ExpansionOrder = Min[ 1, OptionValue[ExpOrder] ];
     
  If[ ! OptionValue[Silent], Print["Expanding lepton vertices..."] ];
  LeptonGaugeVerticesLatex = LeptonGaugeVertices // OutputParametersExpansion;
  LeptonHiggsGaugeVerticesLatex = LeptonHiggsGaugeVertices // OutputParametersExpansion;

  If[ ! OptionValue[Silent], Print["Expanding quark vertices..."] ];
  QuarkGaugeVerticesLatex = QuarkGaugeVertices // OutputParametersExpansion;
  QuarkHiggsGaugeVerticesLatex = QuarkHiggsGaugeVertices // OutputParametersExpansion;
  QuarkGluonVerticesLatex = QuarkGluonVertices // OutputParametersExpansion;

  If[ ! OptionValue[Silent], Print["Expanding gauge and Higgs vertices..."] ];
  GaugeSelfVerticesLatex = GaugeSelfVertices // OutputParametersExpansion;
  GaugeHiggsVerticesLatex = GaugeHiggsVertices // OutputParametersExpansion;
  GluonSelfVerticesLatex = GluonSelfVertices // OutputParametersExpansion;
  GluonHiggsVerticesLatex = GluonHiggsVertices // OutputParametersExpansion;

  If [ SMEFT$RxiGaugeStatus,
       If[ ! OptionValue[Silent], Print["Expanding ghost vertices..."] ];
       GhostVerticesLatex = GhostVertices // OutputParametersExpansion
     ];

  If[ ! OptionValue[Silent], Print["Expanding 4-fermion vertices..."] ];
  FourLeptonVerticesLatex = FourLeptonVertices // OutputParametersExpansion;
  TwoQuarkTwoLeptonVerticesLatex = TwoQuarkTwoLeptonVertices // OutputParametersExpansion;
  FourQuarkVerticesLatex = FourQuarkVertices // OutputParametersExpansion;

  If[ ! OptionValue[Silent], Print["Expanding BL-violating vertices..."] ];
  BLViolatingVerticesLatex = BLViolatingVertices // OutputParametersExpansion;
,
  SMEFT$ExpansionOrder = OptionValue[ExpOrder];

  If[ ! OptionValue[Silent], Print["Expanding lepton vertices..."] ];
  LeptonGaugeVerticesExp = LeptonGaugeVertices // OutputParametersExpansion;
  LeptonHiggsGaugeVerticesExp = LeptonHiggsGaugeVertices // OutputParametersExpansion;

  If[ ! OptionValue[Silent], Print["Expanding quark vertices..."] ];
  QuarkGaugeVerticesExp = QuarkGaugeVertices // OutputParametersExpansion;
  QuarkHiggsGaugeVerticesExp = QuarkHiggsGaugeVertices // OutputParametersExpansion;
  QuarkGluonVerticesExp = QuarkGluonVertices // OutputParametersExpansion;

  If[ ! OptionValue[Silent], Print["Expanding gauge and Higgs vertices..."] ];
  GaugeSelfVerticesExp = GaugeSelfVertices // OutputParametersExpansion;
  GaugeHiggsVerticesExp = GaugeHiggsVertices // OutputParametersExpansion;
  GluonSelfVerticesExp = GluonSelfVertices // OutputParametersExpansion;
  GluonHiggsVerticesExp = GluonHiggsVertices // OutputParametersExpansion;

  If [ SMEFT$RxiGaugeStatus,
       If[ ! OptionValue[Silent], Print["Expanding ghost vertices..."] ];
       GhostVerticesExp = GhostVertices // OutputParametersExpansion
     ];

  If[ ! OptionValue[Silent], Print["Expanding 4-fermion vertices..."] ];
  FourLeptonVerticesExp = FourLeptonVertices // OutputParametersExpansion;
  TwoQuarkTwoLeptonVerticesExp = TwoQuarkTwoLeptonVertices // OutputParametersExpansion;
  FourQuarkVerticesExp = FourQuarkVertices // OutputParametersExpansion;

  If[ ! OptionValue[Silent], Print["Expanding BL-violating vertices..."] ];
  BLViolatingVerticesExp = BLViolatingVertices // OutputParametersExpansion;
];

If[ ! OptionValue[Silent], Print[Style["All vertices expanded, time = ", Bold], TimeUsed[] - CPUTime] ];
       
(* restore general expansion order *)       
SMEFT$ExpansionOrder = order;

(* end of SMEFTExpandVertices *)
]




ReplaceFlavorRotations = Function[x,
(* replace flavor rotation matrices with CKM and PMNS *)
Module[{tmp,ind1,ind2,ff,ff1,ff2,ff3,ff4,a,aa,b,c},

tmp = x;
	
(* apply unitarity of fermion rotations *)
tmp = tmp //. VVL[aa_,b_] Conjugate[VVL[aa_,c_]] -> IndexDelta[b,c];
tmp = tmp //. VLL[aa_,b_] Conjugate[VLL[aa_,c_]] -> IndexDelta[b,c];
tmp = tmp //. VLR[aa_,b_] Conjugate[VLR[aa_,c_]] -> IndexDelta[b,c];
tmp = tmp //. VUL[aa_,b_] Conjugate[VUL[aa_,c_]] -> IndexDelta[b,c];
tmp = tmp //. VUR[aa_,b_] Conjugate[VUR[aa_,c_]] -> IndexDelta[b,c];
tmp = tmp //. VDL[aa_,b_] Conjugate[VDL[aa_,c_]] -> IndexDelta[b,c];
tmp = tmp //. VDR[aa_,b_] Conjugate[VDR[aa_,c_]] -> IndexDelta[b,c];

tmp = tmp //. VVL[b_,aa_] Conjugate[VVL[c_,aa_]] -> IndexDelta[b,c];
tmp = tmp //. VLL[b_,aa_] Conjugate[VLL[c_,aa_]] -> IndexDelta[b,c];
tmp = tmp //. VLR[b_,aa_] Conjugate[VLR[c_,aa_]] -> IndexDelta[b,c];
tmp = tmp //. VUL[b_,aa_] Conjugate[VUL[c_,aa_]] -> IndexDelta[b,c];
tmp = tmp //. VUR[b_,aa_] Conjugate[VUR[c_,aa_]] -> IndexDelta[b,c];
tmp = tmp //. VDL[b_,aa_] Conjugate[VDL[c_,aa_]] -> IndexDelta[b,c];
tmp = tmp //. VDR[b_,aa_] Conjugate[VDR[c_,aa_]] -> IndexDelta[b,c];

(* K=VUL^+ VDL, U = VLL^+ VVL + dim-6 corrections *)
ind1 = Index[Generation,Unique[ff]];
tmp = tmp //. Conjugate[VLL[aa_,b_]] VVL[aa_,c_] -> Ul[b,c] -
              SMEFT$Unitary Lam vev^2 Ul[ind1,c] MBTensor2List["phil3",b,ind1];
ind1 = Index[Generation,Unique[ff]];
tmp = tmp //. Conjugate[VVL[aa_,b_]] VLL[aa_,c_] -> Conjugate[Ul[c,b]] -
              SMEFT$Unitary Lam vev^2 Conjugate[Ul[ind1,b]] Conjugate[ MBTensor2List["phil3",c,ind1] ];
ind1 = Index[Generation,Unique[ff]];
tmp = tmp //. Conjugate[VUL[aa_,b_]] VDL[aa_,c_] -> Kq[b,c] -
              SMEFT$Unitary Lam vev^2 Kq[b,ind1] MBTensor2List["phiq3",ind1,c];
ind1 = Index[Generation,Unique[ff]];
tmp = tmp //. Conjugate[VDL[aa_,b_]] VUL[aa_,c_] -> Conjugate[Kq[c,b]] -
              SMEFT$Unitary Lam vev^2 Conjugate[Kq[c,ind1]] Conjugate[ MBTensor2List["phiq3",ind1,b] ];

tmp

]
(* End of ReplaceFlavorRotations *)
];



ApplyFlavorRotations = Function[x,
(* express 2-fermion WCs in rotated (mass) basis *)
Module[{tmp,i,wcw,wcm,ind1,ind2,ff,ff1,ff2},

tmp = x // Expand // FunctionExpand;

For[i=1, i < Length[Tensor2WC] + 1, i++,
  If[ ! Tensor2WC[[i,5]],
    ind1 = Index[Generation,Unique[ff]];
    ind2 = Index[Generation,Unique[ff]];
    wcw = ToExpression[ SMEFT$WB <> ToString[ Tensor2WC[[i,1]] ] ];
    wcm = ToExpression[ SMEFT$MB <> ToString[ Tensor2WC[[i,1]] ] ];
    tmp = tmp /. wcw[ff1_,ff2_] -> Tensor2WC[[i,2]][ff1,ind1] wcm[ind1,ind2] Conjugate[ Tensor2WC[[i,3]][ff2,ind2] ];
  ];
];

tmp // FunctionExpand

]
(* End of ApplyFlavorRotations *)
];



Apply4fermionRotations = Function[x,
(* express 4-fermion WCs in rotated (mass) basis *)
Module[{tmp,i,wcw,wcm,ind1,ind2,ind3,ind4,ff,ff1,ff2,ff3,ff4},

tmp = x // Expand // FunctionExpand;

For[i=1, i < Length[Tensor4WC] + 1, i++,
  If[ ! Tensor4WC[[i,6]],
    ind1 = Index[Generation,Unique[ff]];
    ind2 = Index[Generation,Unique[ff]];
    ind3 = Index[Generation,Unique[ff]];
    ind4 = Index[Generation,Unique[ff]];
    wcw = ToExpression[ SMEFT$WB <> ToString[ Tensor4WC[[i,1]] ] ];
    wcm = ToExpression[ SMEFT$MB <> ToString[ Tensor4WC[[i,1]] ] ];
    tmp = tmp /. wcw[ff1_,ff2_,ff3_,ff4_] -> Tensor4WC[[i,2]][ff1,ind1]  Tensor4WC[[i,3]][ff3,ind3] wcm[ind1,ind2,ind3,ind4] Conjugate[ Tensor4WC[[i,4]][ff2,ind2] ] Conjugate[ Tensor4WC[[i,5]][ff4,ind4] ];
  ];
];

tmp = tmp // FunctionExpand // Expand;

tmp

]
(* End of Apply4fermionRotations *)
];



ApplyBLFlavorRotations = Function[x,
(* express BL-violating WCs in rotated basis *)
Module[{tmp,i,wcw,wcm,ind1,ind2,ind3,ind4,ff,ff1,ff2,ff3,ff4},

tmp = x // Expand // FunctionExpand;

For[i=1, i < Length[Tensor2WC] + 1, i++,
  If[Tensor2WC[[i,5]],
    ind1 = Index[Generation,Unique[ff]];
    ind2 = Index[Generation,Unique[ff]];
    wcw = ToExpression[ SMEFT$WB <> ToString[ Tensor2WC[[i,1]] ] ];
    wcm = ToExpression[ SMEFT$MB <> ToString[ Tensor2WC[[i,1]] ] ];
    tmp = tmp /. wcw[ff1_,ff2_] -> Conjugate[ Tensor2WC[[i,2]][ff1,ind1] ] wcm[ind1,ind2] Conjugate[ Tensor2WC[[i,3]][ff2,ind2] ];
  ];
];

For[i=1, i < Length[Tensor4WC] + 1, i++,
  If[Tensor4WC[[i,6]],
    ind1 = Index[Generation,Unique[ff]];
    ind2 = Index[Generation,Unique[ff]];
    ind3 = Index[Generation,Unique[ff]];
    ind4 = Index[Generation,Unique[ff]];
    wcw = ToExpression[ SMEFT$WB <> ToString[ Tensor4WC[[i,1]] ] ];
    wcm = ToExpression[ SMEFT$MB <> ToString[ Tensor4WC[[i,1]] ] ];
    tmp = tmp /. wcw[ff1_,ff2_,ff3_,ff4_] -> Conjugate[ Tensor4WC[[i,2]][ff1,ind1] ] Conjugate[ Tensor4WC[[i,3]][ff3,ind3] ] wcm[ind1,ind2,ind3,ind4] Conjugate[ Tensor4WC[[i,4]][ff2,ind2] ] Conjugate[ Tensor4WC[[i,5]][ff4,ind4] ];
  ];
];

tmp = tmp // FunctionExpand // Expand;

tmp

]
(* End of ApplyBLfermionRotations *)
];



Tensor2WCSymmetrize = Function[{WC,cat},
(* symmetrization of 3x3 WC matrices *)
Module[{tmp},

tmp = WC;

If[cat == 2,(* hermitian 3x3 matrix *)
  tmp[[2,1]] = Conjugate[tmp[[1,2]]];
  tmp[[3,1]] = Conjugate[tmp[[1,3]]];
  tmp[[3,2]] = Conjugate[tmp[[2,3]]];
];

If[cat == 3,(* symmetric 3x3 matrix *)
  tmp[[2,1]] = tmp[[1,2]];
  tmp[[3,1]] = tmp[[1,3]];
  tmp[[3,2]] = tmp[[2,3]];
];

tmp

]
(* end of Tensor2WCSymmetrize *)
]



Tensor4WCSymmetrize = Function[{WC,cat},
(* symmetrization of 3x3x3x3 WC matrices, according to their classes
defined by cat argument, see smeft_variables.m *)
Module[{WCR,WCI,tmp},

tmp = WC;

If[ cat == 4 || cat == 5 || cat == 6,

WCR = Re[WC];
WCI = Im[WC];

If[cat == 4,(* two identical XX currents *)
(* Real part *)
  WCR[[1,1,2,1]] = WCR[[1,1,1,2]];
  WCR[[1,1,3,1]] = WCR[[1,1,1,3]];
  WCR[[1,1,3,2]] = WCR[[1,1,2,3]];
  WCR[[1,2,1,1]] = WCR[[1,1,1,2]];
  WCR[[1,3,1,1]] = WCR[[1,1,1,3]];
  WCR[[1,3,1,2]] = WCR[[1,2,1,3]];
  WCR[[1,3,2,1]] = WCR[[1,2,3,1]];
  WCR[[2,1,1,1]] = WCR[[1,1,1,2]];
  WCR[[2,1,1,2]] = WCR[[1,2,2,1]];
  WCR[[2,1,1,3]] = WCR[[1,2,3,1]];
  WCR[[2,1,2,1]] = WCR[[1,2,1,2]];
  WCR[[2,1,2,2]] = WCR[[1,2,2,2]];
  WCR[[2,1,2,3]] = WCR[[1,2,3,2]];
  WCR[[2,1,3,1]] = WCR[[1,2,1,3]];
  WCR[[2,1,3,2]] = WCR[[1,2,2,3]];
  WCR[[2,1,3,3]] = WCR[[1,2,3,3]];
  WCR[[2,2,1,1]] = WCR[[1,1,2,2]];
  WCR[[2,2,1,2]] = WCR[[1,2,2,2]];
  WCR[[2,2,1,3]] = WCR[[1,3,2,2]];
  WCR[[2,2,2,1]] = WCR[[1,2,2,2]];
  WCR[[2,2,3,1]] = WCR[[1,3,2,2]];
  WCR[[2,2,3,2]] = WCR[[2,2,2,3]];
  WCR[[2,3,1,1]] = WCR[[1,1,2,3]];
  WCR[[2,3,1,2]] = WCR[[1,2,2,3]];
  WCR[[2,3,1,3]] = WCR[[1,3,2,3]];
  WCR[[2,3,2,1]] = WCR[[1,2,3,2]];
  WCR[[2,3,2,2]] = WCR[[2,2,2,3]];
  WCR[[2,3,3,1]] = WCR[[1,3,3,2]];
  WCR[[3,1,1,1]] = WCR[[1,1,1,3]];
  WCR[[3,1,1,2]] = WCR[[1,2,3,1]];
  WCR[[3,1,1,3]] = WCR[[1,3,3,1]];
  WCR[[3,1,2,1]] = WCR[[1,2,1,3]];
  WCR[[3,1,2,2]] = WCR[[1,3,2,2]];
  WCR[[3,1,2,3]] = WCR[[1,3,3,2]];
  WCR[[3,1,3,1]] = WCR[[1,3,1,3]];
  WCR[[3,1,3,2]] = WCR[[1,3,2,3]];
  WCR[[3,1,3,3]] = WCR[[1,3,3,3]];
  WCR[[3,2,1,1]] = WCR[[1,1,2,3]];
  WCR[[3,2,1,2]] = WCR[[1,2,3,2]];
  WCR[[3,2,1,3]] = WCR[[1,3,3,2]];
  WCR[[3,2,2,1]] = WCR[[1,2,2,3]];
  WCR[[3,2,2,2]] = WCR[[2,2,2,3]];
  WCR[[3,2,2,3]] = WCR[[2,3,3,2]];
  WCR[[3,2,3,1]] = WCR[[1,3,2,3]];
  WCR[[3,2,3,2]] = WCR[[2,3,2,3]];
  WCR[[3,2,3,3]] = WCR[[2,3,3,3]];
  WCR[[3,3,1,1]] = WCR[[1,1,3,3]];
  WCR[[3,3,1,2]] = WCR[[1,2,3,3]];
  WCR[[3,3,1,3]] = WCR[[1,3,3,3]];
  WCR[[3,3,2,1]] = WCR[[1,2,3,3]];
  WCR[[3,3,2,2]] = WCR[[2,2,3,3]];
  WCR[[3,3,2,3]] = WCR[[2,3,3,3]];
  WCR[[3,3,3,1]] = WCR[[1,3,3,3]];
  WCR[[3,3,3,2]] = WCR[[2,3,3,3]];
(* Imaginary part *)
  WCI[[1,1,1,1]] = 0;
  WCI[[1,1,2,1]] = - WCI[[1,1,1,2]];
  WCI[[1,1,2,2]] = 0;
  WCI[[1,1,3,1]] = - WCI[[1,1,1,3]];
  WCI[[1,1,3,2]] = - WCI[[1,1,2,3]];
  WCI[[1,1,3,3]] = 0;
  WCI[[1,2,1,1]] = WCI[[1,1,1,2]];
  WCI[[1,2,2,1]] = 0;
  WCI[[1,3,1,1]] = WCI[[1,1,1,3]];
  WCI[[1,3,1,2]] = WCI[[1,2,1,3]];
  WCI[[1,3,2,1]] = - WCI[[1,2,3,1]];
  WCI[[1,3,3,1]] = 0;
  WCI[[2,1,1,1]] = - WCI[[1,1,1,2]];
  WCI[[2,1,1,2]] = 0;
  WCI[[2,1,1,3]] = - WCI[[1,2,3,1]];
  WCI[[2,1,2,1]] = - WCI[[1,2,1,2]];
  WCI[[2,1,2,2]] = - WCI[[1,2,2,2]];
  WCI[[2,1,2,3]] = - WCI[[1,2,3,2]];
  WCI[[2,1,3,1]] = - WCI[[1,2,1,3]];
  WCI[[2,1,3,2]] = - WCI[[1,2,2,3]];
  WCI[[2,1,3,3]] = - WCI[[1,2,3,3]];
  WCI[[2,2,1,1]] = 0;
  WCI[[2,2,1,2]] = WCI[[1,2,2,2]];
  WCI[[2,2,1,3]] = WCI[[1,3,2,2]];
  WCI[[2,2,2,1]] = - WCI[[1,2,2,2]];
  WCI[[2,2,2,2]] = 0;
  WCI[[2,2,3,1]] = - WCI[[1,3,2,2]];
  WCI[[2,2,3,2]] = - WCI[[2,2,2,3]];
  WCI[[2,2,3,3]] = 0;
  WCI[[2,3,1,1]] = WCI[[1,1,2,3]];
  WCI[[2,3,1,2]] = WCI[[1,2,2,3]];
  WCI[[2,3,1,3]] = WCI[[1,3,2,3]];
  WCI[[2,3,2,1]] = - WCI[[1,2,3,2]];
  WCI[[2,3,2,2]] = WCI[[2,2,2,3]];
  WCI[[2,3,3,1]] = - WCI[[1,3,3,2]];
  WCI[[2,3,3,2]] = 0;
  WCI[[3,1,1,1]] = - WCI[[1,1,1,3]];
  WCI[[3,1,1,2]] = WCI[[1,2,3,1]];
  WCI[[3,1,1,3]] = 0;
  WCI[[3,1,2,1]] = - WCI[[1,2,1,3]];
  WCI[[3,1,2,2]] = - WCI[[1,3,2,2]];
  WCI[[3,1,2,3]] = - WCI[[1,3,3,2]];
  WCI[[3,1,3,1]] = - WCI[[1,3,1,3]];
  WCI[[3,1,3,2]] = - WCI[[1,3,2,3]];
  WCI[[3,1,3,3]] = - WCI[[1,3,3,3]];
  WCI[[3,2,1,1]] = - WCI[[1,1,2,3]];
  WCI[[3,2,1,2]] = WCI[[1,2,3,2]];
  WCI[[3,2,1,3]] = WCI[[1,3,3,2]];
  WCI[[3,2,2,1]] = - WCI[[1,2,2,3]];
  WCI[[3,2,2,2]] = - WCI[[2,2,2,3]];
  WCI[[3,2,2,3]] = 0;
  WCI[[3,2,3,1]] = - WCI[[1,3,2,3]];
  WCI[[3,2,3,2]] = - WCI[[2,3,2,3]];
  WCI[[3,2,3,3]] = - WCI[[2,3,3,3]];
  WCI[[3,3,1,1]] = 0;
  WCI[[3,3,1,2]] = WCI[[1,2,3,3]];
  WCI[[3,3,1,3]] = WCI[[1,3,3,3]];
  WCI[[3,3,2,1]] = - WCI[[1,2,3,3]];
  WCI[[3,3,2,2]] = 0;
  WCI[[3,3,2,3]] = WCI[[2,3,3,3]];
  WCI[[3,3,3,1]] = - WCI[[1,3,3,3]];
  WCI[[3,3,3,2]] = - WCI[[2,3,3,3]];
  WCI[[3,3,3,3]] = 0;
];

If[cat == 5,(* two independent XX currents *)
(* Real part *)
  WCR[[1,1,2,1]] = WCR[[1,1,1,2]];
  WCR[[1,1,3,1]] = WCR[[1,1,1,3]];
  WCR[[1,1,3,2]] = WCR[[1,1,2,3]];
  WCR[[2,1,1,1]] = WCR[[1,2,1,1]];
  WCR[[2,1,1,2]] = WCR[[1,2,2,1]];
  WCR[[2,1,1,3]] = WCR[[1,2,3,1]];
  WCR[[2,1,2,1]] = WCR[[1,2,1,2]];
  WCR[[2,1,2,2]] = WCR[[1,2,2,2]];
  WCR[[2,1,2,3]] = WCR[[1,2,3,2]];
  WCR[[2,1,3,1]] = WCR[[1,2,1,3]];
  WCR[[2,1,3,2]] = WCR[[1,2,2,3]];
  WCR[[2,1,3,3]] = WCR[[1,2,3,3]];
  WCR[[2,2,2,1]] = WCR[[2,2,1,2]];
  WCR[[2,2,3,1]] = WCR[[2,2,1,3]];
  WCR[[2,2,3,2]] = WCR[[2,2,2,3]];
  WCR[[3,1,1,1]] = WCR[[1,3,1,1]];
  WCR[[3,1,1,2]] = WCR[[1,3,2,1]];
  WCR[[3,1,1,3]] = WCR[[1,3,3,1]];
  WCR[[3,1,2,1]] = WCR[[1,3,1,2]];
  WCR[[3,1,2,2]] = WCR[[1,3,2,2]];
  WCR[[3,1,2,3]] = WCR[[1,3,3,2]];
  WCR[[3,1,3,1]] = WCR[[1,3,1,3]];
  WCR[[3,1,3,2]] = WCR[[1,3,2,3]];
  WCR[[3,1,3,3]] = WCR[[1,3,3,3]];
  WCR[[3,2,1,1]] = WCR[[2,3,1,1]];
  WCR[[3,2,1,2]] = WCR[[2,3,2,1]];
  WCR[[3,2,1,3]] = WCR[[2,3,3,1]];
  WCR[[3,2,2,1]] = WCR[[2,3,1,2]];
  WCR[[3,2,2,2]] = WCR[[2,3,2,2]];
  WCR[[3,2,2,3]] = WCR[[2,3,3,2]];
  WCR[[3,2,3,1]] = WCR[[2,3,1,3]];
  WCR[[3,2,3,2]] = WCR[[2,3,2,3]];
  WCR[[3,2,3,3]] = WCR[[2,3,3,3]];
  WCR[[3,3,2,1]] = WCR[[3,3,1,2]];
  WCR[[3,3,3,1]] = WCR[[3,3,1,3]];
  WCR[[3,3,3,2]] = WCR[[3,3,2,3]];
(* Imaginary part *)
  WCI[[1,1,1,1]] = 0;
  WCI[[1,1,2,1]] = - WCI[[1,1,1,2]];
  WCI[[1,1,2,2]] = 0;
  WCI[[1,1,3,1]] = - WCI[[1,1,1,3]];
  WCI[[1,1,3,2]] = - WCI[[1,1,2,3]];
  WCI[[1,1,3,3]] = 0;
  WCI[[2,1,1,1]] = - WCI[[1,2,1,1]];
  WCI[[2,1,1,2]] = - WCI[[1,2,2,1]];
  WCI[[2,1,1,3]] = - WCI[[1,2,3,1]];
  WCI[[2,1,2,1]] = - WCI[[1,2,1,2]];
  WCI[[2,1,2,2]] = - WCI[[1,2,2,2]];
  WCI[[2,1,2,3]] = - WCI[[1,2,3,2]];
  WCI[[2,1,3,1]] = - WCI[[1,2,1,3]];
  WCI[[2,1,3,2]] = - WCI[[1,2,2,3]];
  WCI[[2,1,3,3]] = - WCI[[1,2,3,3]];
  WCI[[2,2,1,1]] = 0;
  WCI[[2,2,2,1]] = - WCI[[2,2,1,2]];
  WCI[[2,2,2,2]] = 0;
  WCI[[2,2,3,1]] = - WCI[[2,2,1,3]];
  WCI[[2,2,3,2]] = - WCI[[2,2,2,3]];
  WCI[[2,2,3,3]] = 0;
  WCI[[3,1,1,1]] = - WCI[[1,3,1,1]];
  WCI[[3,1,1,2]] = - WCI[[1,3,2,1]];
  WCI[[3,1,1,3]] = - WCI[[1,3,3,1]];
  WCI[[3,1,2,1]] = - WCI[[1,3,1,2]];
  WCI[[3,1,2,2]] = - WCI[[1,3,2,2]];
  WCI[[3,1,2,3]] = - WCI[[1,3,3,2]];
  WCI[[3,1,3,1]] = - WCI[[1,3,1,3]];
  WCI[[3,1,3,2]] = - WCI[[1,3,2,3]];
  WCI[[3,1,3,3]] = - WCI[[1,3,3,3]];
  WCI[[3,2,1,1]] = - WCI[[2,3,1,1]];
  WCI[[3,2,1,2]] = - WCI[[2,3,2,1]];
  WCI[[3,2,1,3]] = - WCI[[2,3,3,1]];
  WCI[[3,2,2,1]] = - WCI[[2,3,1,2]];
  WCI[[3,2,2,2]] = - WCI[[2,3,2,2]];
  WCI[[3,2,2,3]] = - WCI[[2,3,3,2]];
  WCI[[3,2,3,1]] = - WCI[[2,3,1,3]];
  WCI[[3,2,3,2]] = - WCI[[2,3,2,3]];
  WCI[[3,2,3,3]] = - WCI[[2,3,3,3]];
  WCI[[3,3,1,1]] = 0;
  WCI[[3,3,2,1]] = - WCI[[3,3,1,2]];
  WCI[[3,3,2,2]] = 0;
  WCI[[3,3,3,1]] = - WCI[[3,3,1,3]];
  WCI[[3,3,3,2]] = - WCI[[3,3,2,3]];
  WCI[[3,3,3,3]] = 0;
];

If[cat == 6,(* two identical XX currents special case Cee *)
(* Real part *)
  WCR[[1,1,2,1]] = WCR[[1,1,1,2]];
  WCR[[1,1,3,1]] = WCR[[1,1,1,3]];
  WCR[[1,1,3,2]] = WCR[[1,1,2,3]];
  WCR[[1,2,1,1]] = WCR[[1,1,1,2]];
  WCR[[1,2,2,1]] = WCR[[1,1,2,2]];
  WCR[[1,2,3,1]] = WCR[[1,1,2,3]];
  WCR[[1,3,1,1]] = WCR[[1,1,1,3]];
  WCR[[1,3,1,2]] = WCR[[1,2,1,3]];
  WCR[[1,3,2,1]] = WCR[[1,1,2,3]];
  WCR[[1,3,2,2]] = WCR[[1,2,2,3]];
  WCR[[1,3,3,1]] = WCR[[1,1,3,3]];
  WCR[[1,3,3,2]] = WCR[[1,2,3,3]];
  WCR[[2,1,1,1]] = WCR[[1,1,1,2]];
  WCR[[2,1,1,2]] = WCR[[1,1,2,2]];
  WCR[[2,1,1,3]] = WCR[[1,1,2,3]];
  WCR[[2,1,2,1]] = WCR[[1,2,1,2]];
  WCR[[2,1,2,2]] = WCR[[1,2,2,2]];
  WCR[[2,1,2,3]] = WCR[[1,2,3,2]];
  WCR[[2,1,3,1]] = WCR[[1,2,1,3]];
  WCR[[2,1,3,2]] = WCR[[1,2,2,3]];
  WCR[[2,1,3,3]] = WCR[[1,2,3,3]];
  WCR[[2,2,1,1]] = WCR[[1,1,2,2]];
  WCR[[2,2,1,2]] = WCR[[1,2,2,2]];
  WCR[[2,2,1,3]] = WCR[[1,2,2,3]];
  WCR[[2,2,2,1]] = WCR[[1,2,2,2]];
  WCR[[2,2,3,1]] = WCR[[1,2,2,3]];
  WCR[[2,2,3,2]] = WCR[[2,2,2,3]];
  WCR[[2,3,1,1]] = WCR[[1,1,2,3]];
  WCR[[2,3,1,2]] = WCR[[1,2,2,3]];
  WCR[[2,3,1,3]] = WCR[[1,3,2,3]];
  WCR[[2,3,2,1]] = WCR[[1,2,3,2]];
  WCR[[2,3,2,2]] = WCR[[2,2,2,3]];
  WCR[[2,3,3,1]] = WCR[[1,2,3,3]];
  WCR[[2,3,3,2]] = WCR[[2,2,3,3]];
  WCR[[3,1,1,1]] = WCR[[1,1,1,3]];
  WCR[[3,1,1,2]] = WCR[[1,1,2,3]];
  WCR[[3,1,1,3]] = WCR[[1,1,3,3]];
  WCR[[3,1,2,1]] = WCR[[1,2,1,3]];
  WCR[[3,1,2,2]] = WCR[[1,2,2,3]];
  WCR[[3,1,2,3]] = WCR[[1,2,3,3]];
  WCR[[3,1,3,1]] = WCR[[1,3,1,3]];
  WCR[[3,1,3,2]] = WCR[[1,3,2,3]];
  WCR[[3,1,3,3]] = WCR[[1,3,3,3]];
  WCR[[3,2,1,1]] = WCR[[1,1,2,3]];
  WCR[[3,2,1,2]] = WCR[[1,2,3,2]];
  WCR[[3,2,1,3]] = WCR[[1,2,3,3]];
  WCR[[3,2,2,1]] = WCR[[1,2,2,3]];
  WCR[[3,2,2,2]] = WCR[[2,2,2,3]];
  WCR[[3,2,2,3]] = WCR[[2,2,3,3]];
  WCR[[3,2,3,1]] = WCR[[1,3,2,3]];
  WCR[[3,2,3,2]] = WCR[[2,3,2,3]];
  WCR[[3,2,3,3]] = WCR[[2,3,3,3]];
  WCR[[3,3,1,1]] = WCR[[1,1,3,3]];
  WCR[[3,3,1,2]] = WCR[[1,2,3,3]];
  WCR[[3,3,1,3]] = WCR[[1,3,3,3]];
  WCR[[3,3,2,1]] = WCR[[1,2,3,3]];
  WCR[[3,3,2,2]] = WCR[[2,2,3,3]];
  WCR[[3,3,2,3]] = WCR[[2,3,3,3]];
  WCR[[3,3,3,1]] = WCR[[1,3,3,3]];
  WCR[[3,3,3,2]] = WCR[[2,3,3,3]];
(* Imaginary part *)
  WCI[[1,1,1,1]] = 0;
  WCI[[1,1,2,1]] = - WCI[[1,1,1,2]];
  WCI[[1,1,2,2]] = 0;
  WCI[[1,1,3,1]] = - WCI[[1,1,1,3]];
  WCI[[1,1,3,2]] = - WCI[[1,1,2,3]];
  WCI[[1,1,3,3]] = 0;
  WCI[[1,2,1,1]] = WCI[[1,1,1,2]];
  WCI[[1,2,2,1]] = 0;
  WCI[[1,2,3,1]] = - WCI[[1,1,2,3]];
  WCI[[1,3,1,1]] = WCI[[1,1,1,3]];
  WCI[[1,3,1,2]] = WCI[[1,2,1,3]];
  WCI[[1,3,2,1]] = WCI[[1,1,2,3]];
  WCI[[1,3,2,2]] = WCI[[1,2,2,3]];
  WCI[[1,3,3,1]] = 0;
  WCI[[1,3,3,2]] = WCI[[1,2,3,3]];
  WCI[[2,1,1,1]] = - WCI[[1,1,1,2]];
  WCI[[2,1,1,2]] = 0;
  WCI[[2,1,1,3]] = WCI[[1,1,2,3]];
  WCI[[2,1,2,1]] = - WCI[[1,2,1,2]];
  WCI[[2,1,2,2]] = - WCI[[1,2,2,2]];
  WCI[[2,1,2,3]] = - WCI[[1,2,3,2]];
  WCI[[2,1,3,1]] = - WCI[[1,2,1,3]];
  WCI[[2,1,3,2]] = - WCI[[1,2,2,3]];
  WCI[[2,1,3,3]] = - WCI[[1,2,3,3]];
  WCI[[2,2,1,1]] = 0;
  WCI[[2,2,1,2]] = WCI[[1,2,2,2]];
  WCI[[2,2,1,3]] = WCI[[1,2,2,3]];
  WCI[[2,2,2,1]] = - WCI[[1,2,2,2]];
  WCI[[2,2,2,2]] = 0;
  WCI[[2,2,3,1]] = - WCI[[1,2,2,3]];
  WCI[[2,2,3,2]] = - WCI[[2,2,2,3]];
  WCI[[2,2,3,3]] = 0;
  WCI[[2,3,1,1]] = WCI[[1,1,2,3]];
  WCI[[2,3,1,2]] = WCI[[1,2,2,3]];
  WCI[[2,3,1,3]] = WCI[[1,3,2,3]];
  WCI[[2,3,2,1]] = - WCI[[1,2,3,2]];
  WCI[[2,3,2,2]] = WCI[[2,2,2,3]];
  WCI[[2,3,3,1]] = - WCI[[1,2,3,3]];
  WCI[[2,3,3,2]] = 0;
  WCI[[3,1,1,1]] = - WCI[[1,1,1,3]];
  WCI[[3,1,1,2]] = - WCI[[1,1,2,3]];
  WCI[[3,1,1,3]] = 0;
  WCI[[3,1,2,1]] = - WCI[[1,2,1,3]];
  WCI[[3,1,2,2]] = - WCI[[1,2,2,3]];
  WCI[[3,1,2,3]] = - WCI[[1,2,3,3]];
  WCI[[3,1,3,1]] = - WCI[[1,3,1,3]];
  WCI[[3,1,3,2]] = - WCI[[1,3,2,3]];
  WCI[[3,1,3,3]] = - WCI[[1,3,3,3]];
  WCI[[3,2,1,1]] = - WCI[[1,1,2,3]];
  WCI[[3,2,1,2]] = WCI[[1,2,3,2]];
  WCI[[3,2,1,3]] = WCI[[1,2,3,3]];
  WCI[[3,2,2,1]] = - WCI[[1,2,2,3]];
  WCI[[3,2,2,2]] = - WCI[[2,2,2,3]];
  WCI[[3,2,2,3]] = 0;
  WCI[[3,2,3,1]] = - WCI[[1,3,2,3]];
  WCI[[3,2,3,2]] = - WCI[[2,3,2,3]];
  WCI[[3,2,3,3]] = - WCI[[2,3,3,3]];
  WCI[[3,3,1,1]] = 0;
  WCI[[3,3,1,2]] = WCI[[1,2,3,3]];
  WCI[[3,3,1,3]] = WCI[[1,3,3,3]];
  WCI[[3,3,2,1]] = - WCI[[1,2,3,3]];
  WCI[[3,3,2,2]] = 0;
  WCI[[3,3,2,3]] = WCI[[2,3,3,3]];
  WCI[[3,3,3,1]] = - WCI[[1,3,3,3]];
  WCI[[3,3,3,2]] = - WCI[[2,3,3,3]];
  WCI[[3,3,3,3]] = 0;
];

tmp = WCR + I WCI;

];

If[cat == 7,(* B- violating special case qque *)
  tmp[[2,1,1,1]] = tmp[[1,2,1,1]];
  tmp[[2,1,1,2]] = tmp[[1,2,1,2]];
  tmp[[2,1,1,3]] = tmp[[1,2,1,3]];
  tmp[[2,1,2,1]] = tmp[[1,2,2,1]];
  tmp[[2,1,2,2]] = tmp[[1,2,2,2]];
  tmp[[2,1,2,3]] = tmp[[1,2,2,3]];
  tmp[[2,1,3,1]] = tmp[[1,2,3,1]];
  tmp[[2,1,3,2]] = tmp[[1,2,3,2]];
  tmp[[2,1,3,3]] = tmp[[1,2,3,3]];
  tmp[[3,1,1,1]] = tmp[[1,3,1,1]];
  tmp[[3,1,1,2]] = tmp[[1,3,1,2]];
  tmp[[3,1,1,3]] = tmp[[1,3,1,3]];
  tmp[[3,1,2,1]] = tmp[[1,3,2,1]];
  tmp[[3,1,2,2]] = tmp[[1,3,2,2]];
  tmp[[3,1,2,3]] = tmp[[1,3,2,3]];
  tmp[[3,1,3,1]] = tmp[[1,3,3,1]];
  tmp[[3,1,3,2]] = tmp[[1,3,3,2]];
  tmp[[3,1,3,3]] = tmp[[1,3,3,3]];
  tmp[[3,2,1,1]] = tmp[[2,3,1,1]];
  tmp[[3,2,1,2]] = tmp[[2,3,1,2]];
  tmp[[3,2,1,3]] = tmp[[2,3,1,3]];
  tmp[[3,2,2,1]] = tmp[[2,3,2,1]];
  tmp[[3,2,2,2]] = tmp[[2,3,2,2]];
  tmp[[3,2,2,3]] = tmp[[2,3,2,3]];
  tmp[[3,2,3,1]] = tmp[[2,3,3,1]];
  tmp[[3,2,3,2]] = tmp[[2,3,3,2]];
  tmp[[3,2,3,3]] = tmp[[2,3,3,3]];
];

If[cat == 8,(* B- violating special case qqql *)
  tmp[[2,1,1,1]] = tmp[[1,1,2,1]];
  tmp[[2,1,1,2]] = tmp[[1,1,2,2]];
  tmp[[2,1,1,3]] = tmp[[1,1,2,3]];
  tmp[[2,2,1,1]] = tmp[[1,2,2,1]];
  tmp[[2,2,1,2]] = tmp[[1,2,2,2]];
  tmp[[2,2,1,3]] = tmp[[1,2,2,3]];
  tmp[[3,1,1,1]] = tmp[[1,1,3,1]];
  tmp[[3,1,1,2]] = tmp[[1,1,3,2]];
  tmp[[3,1,1,3]] = tmp[[1,1,3,3]];
  tmp[[3,1,2,1]] = tmp[[2,3,1,1]] + tmp[[2,1,3,1]] - tmp[[1,3,2,1]];
  tmp[[3,1,2,2]] = tmp[[2,3,1,2]] + tmp[[2,1,3,2]] - tmp[[1,3,2,2]];
  tmp[[3,1,2,3]] = tmp[[2,3,1,3]] + tmp[[2,1,3,3]] - tmp[[1,3,2,3]];
  tmp[[3,2,1,1]] = tmp[[1,3,2,1]] + tmp[[1,2,3,1]] - tmp[[2,3,1,1]];
  tmp[[3,2,1,2]] = tmp[[1,3,2,2]] + tmp[[1,2,3,2]] - tmp[[2,3,1,2]];
  tmp[[3,2,1,3]] = tmp[[1,3,2,3]] + tmp[[1,2,3,3]] - tmp[[2,3,1,3]];
  tmp[[3,2,2,1]] = tmp[[2,2,3,1]];
  tmp[[3,2,2,2]] = tmp[[2,2,3,2]];
  tmp[[3,2,2,3]] = tmp[[2,2,3,3]];
  tmp[[3,3,1,1]] = tmp[[1,3,3,1]];
  tmp[[3,3,1,2]] = tmp[[1,3,3,2]];
  tmp[[3,3,1,3]] = tmp[[1,3,3,3]];
  tmp[[3,3,2,1]] = tmp[[2,3,3,1]];
  tmp[[3,3,2,2]] = tmp[[2,3,3,2]];
  tmp[[3,3,2,3]] = tmp[[2,3,3,3]];
];

tmp

]
(* end of Tensor4WCSymmetrize *)
]




SMEFTFeynmanRules = Function[{},
(* calculate SMEFT Feynman rules, sector by sector *)

LeptonEWInteractions[SMEFT$MaxParticles];
DeltaLTwoInteractions[SMEFT$MaxParticles];

If[ SMEFT$MajoranaNeutrino,
  LeptonGaugeVertices = SymmetrizeNeutrino2Current[LeptonGaugeVertices];
  LeptonHiggsGaugeVertices = MergeVertices[DeltaLTwoVertices, LeptonHiggsGaugeVertices];
  LeptonHiggsGaugeVertices = SymmetrizeNeutrino2Current[LeptonHiggsGaugeVertices];
];

QuarkEWInteractions[SMEFT$MaxParticles];
QuarkQCDInteractions[SMEFT$MaxParticles];

HiggsEWInteractions[SMEFT$MaxParticles];
HiggsQCDInteractions[SMEFT$MaxParticles];

FourFermionInteractions[SMEFT$MaxParticles];
BLViolating4FermionInteractions[SMEFT$MaxParticles];

If[SMEFT$RxiGaugeStatus, GhostInteractions[] ];

(* end of SMEFTFeynmanRules *)
]





SymmetrizeNeutrino2Current = Function[{x},
(* symmetrize vertices for Majorana neutrinos *)
Module[{a,b,c,swap,i,frule,symlist,tmp},

symlist = {};
tmp = x;
       
For[i=1, i < Length[tmp]+1, i++,
(* find symmetrized form of vlbar vl vertices *)
  If[ SubsetQ[ Flatten[ Take[ tmp[[i,1]], All, 1] ], {vlbar,vl} ], 
    Clear[swap]; 
(* swap external indices *)
    frule = tmp[[i]] /. Ext[i1] -> swap;
    frule = frule /. Ext[i2] -> Ext[i1];
    frule = frule /. swap -> Ext[i2];
(* but swap back external spin indices *)
    frule = frule /. Index[Spin,Ext[i1]] -> Index[Spin,swap];
    frule = frule /. Index[Spin,Ext[i2]] -> Index[Spin,Ext[i1]];
    frule = frule /. Index[Spin,swap]    -> Index[Spin,Ext[i2]];
    frule = frule /. TensDot[SlashedP[a_], ProjM][b_,c_] ->    TensDot[SlashedP[a], swap][b,c];
    frule = frule /. TensDot[SlashedP[a_], ProjP][b_,c_] ->  - TensDot[SlashedP[a], ProjM][b,c];
    frule = frule /. TensDot[SlashedP[a_], swap][b_,c_]  ->  - TensDot[SlashedP[a], ProjP][b,c];
    frule = frule /. TensDot[Ga[a_], ProjM][b_,c_] ->    TensDot[Ga[a], swap][b,c];
    frule = frule /. TensDot[Ga[a_], ProjP][b_,c_] ->  - TensDot[Ga[a], ProjM][b,c];
    frule = frule /. TensDot[Ga[a_], swap][b_,c_]  ->  - TensDot[Ga[a], ProjP][b,c];
    AppendTo[symlist,frule];
  ];
(*  change vl vl and vlbar vlbar vertices into vlbar vl, they are identical for Majorana neutrinos *)
  If[ Flatten[ Take[ tmp[[i,1]], All, 1] ][[1;;2]] === {vl,vl}, tmp[[i,1,1,1]] = vlbar ];
  If[ Flatten[ Take[ tmp[[i,1]], All, 1] ][[1;;2]] === {vlbar,vlbar}, tmp[[i,1,2,1]] = vl ];
(*  change lbar vlbar vertex into lbar vl, identical for Majorana neutrinos *)
  If[ Flatten[ Take[ tmp[[i,1]], All, 1] ][[1;;2]] === {lbar,vlbar}, tmp[[i,1,2,1]] = vl ];
(*  change l vl vertex into vlbar l, identical for Majorana neutrinos *)
  If[ Flatten[ Take[ tmp[[i,1]], All, 1] ][[1;;2]] === {l,vl}, tmp[[i,1,1,1]] = vlbar; tmp[[i,1,2,1]] = l ];

];

MergeVertices[tmp,symlist] // Simplify

]
(* end of SymmetrizeNeutrino2Current *)
]


