(* SmeftFR v3.0 package *)
(* Calculation of Feynman rules for gauge and Higgs couplings in
physical basis *)

TwoFermionFlavorRotate = Function[{x},
(* 2 fermion current flavor rotation and expansion *)
ReplaceFlavorRotations[x // Expand] // SMEFTExpReduce
(* end of TwoFermionFlavorRotate *)
]



LeptonEWInteractions = Function[{MaxLeg},
(* 2-lepton EW couplings *)
Module[{aa, bb, tmp, LGG, LGH},

(* leptonic vertices only *)

tmp = SMEFT$LGferm + SMEFT$LGferm6/. G[___] -> 0 /. uq[___] -> 0 /. dq[___] -> 0;
(* express Yukawa's in term of fermion masses *)
tmp = tmp /. yl[aa_,bb_] -> SMEFT$YL[aa,bb] // FunctionExpand;
(* pure gauge vertices *)
LGG = tmp /. H -> 0 /. G0->0 /. GP->0 /. GPbar->0;
(* gauge-Higgs part *)
LGH = tmp - LGG;

LeptonGaugeLagrangian = TwoFermionFlavorRotate[LGG] /. SMEFT$Identities // OptimizeIndex;
LeptonHiggsGaugeLagrangian = TwoFermionFlavorRotate[LGH] /. SMEFT$Identities // OptimizeIndex;
LeptonHiggsGaugeLagrangian = LeptonHiggsGaugeLagrangian /. Conjugate[fml[aa_, bb_]] -> fml[bb, aa] ;

Print[Style["Calculating leptonic EW vertices...",Bold]];
LeptonGaugeVertices = FeynmanRules[LeptonGaugeLagrangian , MaxParticles -> MaxLeg];
LeptonHiggsGaugeVertices = FeynmanRules[LeptonHiggsGaugeLagrangian, MaxParticles -> MaxLeg];

(* LeptonHiggsGaugeVertices = LeptonHiggsGaugeVertices /. fml[a_, b_]  -> mldiag[a] IndexDelta[a,b] // Simplify;*)

];
(* end of LeptonEWInteractions*)
];




QuarkEWInteractions = Function[{MaxLeg},
(* 2-quark EW couplings *)
Module[{ aa, bb, tmp, LGG, LGH},

(* quark vertices only *)
tmp = SMEFT$LGferm /. G[___] -> 0 /. l[___] -> 0 /. vl[___] -> 0;
(* express Yukawa's in term of fermion masses *)
tmp = tmp /. yu[aa_,bb_] -> SMEFT$YU[aa,bb];
tmp = tmp /. yd[aa_,bb_] -> SMEFT$YD[aa,bb] // FunctionExpand;
(* pure gauge vertices *)
LGG = tmp /. H -> 0 /. G0->0 /. GP->0 /. GPbar->0;
(* gauge-Higgs part *)
LGH = tmp - LGG;

QuarkGaugeLagrangian      = TwoFermionFlavorRotate[LGG]  /. SMEFT$Identities // OptimizeIndex;
QuarkHiggsGaugeLagrangian = TwoFermionFlavorRotate[LGH]  /. SMEFT$Identities // OptimizeIndex;
QuarkHiggsGaugeLagrangian = QuarkHiggsGaugeLagrangian /. Conjugate[fmd[aa_, bb_]] -> fmd[bb, aa] /.
                                                         Conjugate[fmu[aa_, bb_]] -> fmu[bb, aa];

Print[Style["Calculating quark EW vertices...",Bold]];
QuarkGaugeVertices = FeynmanRules[QuarkGaugeLagrangian, MaxParticles -> MaxLeg];
QuarkHiggsGaugeVertices = FeynmanRules[QuarkHiggsGaugeLagrangian, MaxParticles -> MaxLeg];

(* QuarkHiggsGaugeVertices = QuarkHiggsGaugeVertices /. fmd[a_, b_]  -> mddiag[a] IndexDelta[a,b] /.
                                                     fmu[a_, b_]  -> mudiag[a] IndexDelta[a,b] // Simplify; *)

];
(* end of QuarkEWInteractions*)
];




QuarkQCDInteractions = Function[{MaxLeg},
(* 2-fermion QCD couplings *)
Module[{LGlF},

LGlF = SMEFT$LGferm - (SMEFT$LGferm /. G[___]->0) // FunctionExpand;

QuarkGluonLagrangian = TwoFermionFlavorRotate[LGlF] /. SMEFT$Identities // Expand;
QuarkGluonLagrangian = OptimizeIndex[QuarkGluonLagrangian] // SMEFTExpReduce;

Print[Style["Calculating quark QCD vertices...",Bold]];
QuarkGluonVertices =  FeynmanRules[QuarkGluonLagrangian, MaxParticles -> MaxLeg];

];
(* end of QuarkQCDInteractions*)
];



DeltaLTwoInteractions = Function[{MaxLeg},
(*  Delta L=2 violating couplings *)
Block[{aa,bb,vev},
				       
Print[Style["Calculating Delta L=2 violating vertices...",Bold]];

SetAttributes[fmv, Orderless];

DeltaLTwoLagrangian = If[ MemberQ[ SMEFT$OperatorList, "vv" ], LQvv, 0] // Expand;
DeltaLTwoLagrangian = DeltaLTwoLagrangian /. ToExpression[ SMEFT$MB <>"vv" ][aa_, bb_] -> - fmv[aa,bb]/vev^2;
DeltaLTwoLagrangian = ReplaceFlavorRotations[DeltaLTwoLagrangian]  /. SMEFT$Identities // OptimizeIndex;
DeltaLTwoVertices = FeynmanRules[DeltaLTwoLagrangian, MaxParticles -> MaxLeg];

];      
(* end of DeltaLTwoInteractions *)
];



HiggsEWInteractions = Function[{MaxLeg},
(* EW gauge and Higgs couplings *)
Module[{a, b, c, A, B, C, k, l, m, LGGH, ind1, ind2, ind3, lgind,
epstmp, lggh },

LGGH = SMEFT$LGeff + SMEFT$LGeff6 /. G[___] -> 0  /. SMEFT$Identities // OptimizeIndex;
LGGH = LGGH // SMEFTExpReduce;

GaugeSelfLagrangian = LGGH /. H->0 /. GP->0 /. GPbar->0 /. G0->0;
GaugeHiggsLagrangian = LGGH - GaugeSelfLagrangian // Expand;

Print[Style["Calculating EW gauge-gauge vertices...",Bold]];
GaugeSelfVertices = FeynmanRules[GaugeSelfLagrangian, MaxParticles -> MaxLeg];
Print[Style["Calculating EW gauge+Higgs vertices...",Bold]];
GaugeHiggsVertices = FeynmanRules[GaugeHiggsLagrangian, MaxParticles -> MaxLeg];

(* some manual index reordering/simplification *)
lggh = {GaugeSelfVertices,GaugeHiggsVertices};

ind1 = Unique[lgind];
ind2 = Unique[lgind];
ind3 = Unique[lgind];

lggh = lggh /. Eps -> epstmp;

lggh = lggh /. epstmp[Index[Lorentz,Ext[k_]], a_, b_, c_] FV[A_, a_] FV[B_, b_] FV[C_, c_] ->
  epstmp[Index[Lorentz,Ext[k]],ind1,ind2,ind3] FV[A,ind1] FV[B,ind2] FV[C,ind3];

lggh = lggh /. epstmp[Index[Lorentz,Ext[k_]], Index[Lorentz,Ext[l_]], a_, b_] FV[A_, a_] FV[B_, b_] ->
  epstmp[Index[Lorentz,Ext[k]], Index[Lorentz,Ext[l]],ind1,ind2] FV[A,ind1] FV[B,ind2];

lggh = lggh /. epstmp[Index[Lorentz,Ext[k_]], Index[Lorentz,Ext[l_]], Index[Lorentz,Ext[m_]], a_] FV[A_, a_] ->
  epstmp[Index[Lorentz,Ext[k]], Index[Lorentz,Ext[l]], Index[Lorentz,Ext[m]],ind1] FV[A,ind1];

lggh = lggh /. ind1-> Index[Lorentz,\[Alpha]1] /. ind2-> Index[Lorentz,\[Beta]1] /. ind3-> Index[Lorentz,\[Gamma]1];

lggh = lggh /. epstmp -> Eps;

GaugeSelfVertices = lggh[[1]];
GaugeHiggsVertices = lggh[[2]] // Simplify;

];
(* end of HiggsEWInteractions *)
];




HiggsQCDInteractions = Function[{MaxLeg},
(* QCD gauge and Higgs couplings *)
Module[{a, b, c, A, B, C, k, l, m, tmp, LGall, LGluon, lgind, ind1,
ind2, ind3, epstmp},

LGall = SMEFT$LGeff + SMEFT$LGeff6 // FunctionExpand;
LGall = LGall  /. SMEFT$Identities // OptimizeIndex;
tmp = LGall /. G[___]->0;
LGluon = Expand[LGall - tmp] // SMEFTExpReduce;

GluonSelfLagrangian = LGluon /. H->0 /. GP->0 /. GPbar->0 /. G0->0;
GluonHiggsLagrangian = LGluon - GluonSelfLagrangian // Expand;

Print[Style["Calculating QCD gauge+Higgs vertices...",Bold]];
GluonSelfVertices  = FeynmanRules[GluonSelfLagrangian, MaxParticles -> MaxLeg];
GluonHiggsVertices = FeynmanRules[GluonHiggsLagrangian, MaxParticles -> MaxLeg];

(* manual index simplifications *)
tmp = {GluonSelfVertices,GluonHiggsVertices};

ind1 = Unique[lgind];
ind2 = Unique[lgind];
ind3 = Unique[lgind];

tmp = tmp /. Eps -> epstmp;

tmp = tmp /. epstmp[Index[Lorentz,Ext[k_]], a_, b_, c_] FV[A_, a_] FV[B_, b_] FV[C_, c_] ->
  epstmp[Index[Lorentz,Ext[k]],ind1,ind2,ind3] FV[A,ind1] FV[B,ind2] FV[C,ind3];

tmp = tmp /. epstmp[Index[Lorentz,Ext[k_]], Index[Lorentz,Ext[l_]], a_, b_] FV[A_, a_] FV[B_, b_] ->
  epstmp[Index[Lorentz,Ext[k]], Index[Lorentz,Ext[l]],ind1,ind2] FV[A,ind1] FV[B,ind2];

tmp = tmp /. epstmp[Index[Lorentz,Ext[k_]], Index[Lorentz,Ext[l_]], Index[Lorentz,Ext[m_]], a_] FV[A_, a_] ->
  epstmp[Index[Lorentz,Ext[k]], Index[Lorentz,Ext[l]], Index[Lorentz,Ext[m]],ind1] FV[A,ind1];

tmp = tmp /. ind1-> Index[Lorentz,\[Alpha]1] /. ind2-> Index[Lorentz,\[Beta]1] /. ind3-> Index[Lorentz,\[Gamma]1];

tmp = tmp /. epstmp -> Eps;

GluonSelfVertices = tmp[[1]];
GluonHiggsVertices = tmp[[2]] // Simplify;

];
(* end of HiggsQCDInteractions *)
];


