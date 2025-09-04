(* SmeftFR v3.0 package *)
(* rediagonalization of SMEFT Lagrangian, field and couplings redefinitions *)


SMEFTGaugeNormalization = Function[{},
(* gauge coupling normalization  g -> gnorm G; add higher dimension operators if necessary *)

SMEFT$g1norm = Sqrt[ 1 - 2 Lam vev^2 MBScalarList["phiB"] - Lam^2 vev^4 MBScalarList["B2phi4n1"] ] // SMEFTExpOrder;

SMEFT$gwnorm = Sqrt[ 1 - 2 Lam vev^2 MBScalarList["phiW"] - Lam^2 vev^4 MBScalarList["W2phi4n1"] ] // SMEFTExpOrder;

SMEFT$gsnorm = Sqrt[ 1 - 2 Lam vev^2 MBScalarList["phiG"] - Lam^2 vev^4 MBScalarList["G2phi4n1"] ] // SMEFTExpOrder;

SMEFT$G1 = G1;
SMEFT$GW = GW;
SMEFT$GS = GS;

(* end of SMEFTGaugeNormalization *)
]


SMEFTFindMassBasis = Function[{},
(* find transformations to mass basis, sector by sector *)

SMEFTGaugeNormalization[];

HiggsSectorDiagonalization[];
GaugeSectorDiagonalization[];
FermionSectorDiagonalization[];

UpdateSpecialParameters[];
UpdateNormalizationConstants[]; 
			      
(* end of SMEFTFindMassBasis *)
]




HiggsSectorDiagonalization = Function[{},
(* HIGGS SECTOR REDIAGONALIZATION *)
Block[{vev},Module[{aa, c1, c2, LHeff, KinH, tmp, VGHiggs, VHiggs, dV, dV0, dV1, vev0,
delv1, delv2, dvev, MGH2, MG0, MGP, VEV},

(* one of standard parameters *)
SMEFT$hlambda = hlambda;

(* effective pure Higgs Lagrangian, Lam is actually 1/Lambda^2! *)
LHeff = SMEFT$LGeff /. A[___] -> 0 /. Z[___] -> 0 /. W[___] -> 0 /. Wbar[___] -> 0 /. G[___] ->0 // Expand;

(* effective Higgs potential and derivative *)
VGHiggs = - LHeff /. del[H_, ___]->0;
VHiggs = VGHiggs /. G0->0 /. GP->0 /. GPbar -> 0;
dV = D[VHiggs,H];

(* we want dV vanish for H=0, so H=0 is true vaccum; this gives
equation for vev *)
dV0 = dV/vev /. H->0 /. Hnorm->1 // Simplify;
(* solve for 0th order in Lambda, choose positive non-0 solution *)
tmp = vev /. Solve[ (dV0 /. Lam -> 0) == 0, vev ];
vev0 = tmp[[2]];
(* higher order Lambda correction *)
dV1 = (dV0 /. vev -> vev0 + Lam delv1 + Lam^2 delv2) // SMEFTExpOrder;
If [SMEFT$ExpansionOrder == 0, dvev = 0];
If [SMEFT$ExpansionOrder == 1, dvev = Lam delv1 /. Solve[ dV1 == 0, delv1 ][[1]] ];
If [SMEFT$ExpansionOrder == 2,
  c1 = Coefficient[dV1, Lam];
  c2 = Coefficient[dV1, Lam^2];
  dvev = Lam delv1 + Lam^2 delv2 /. Solve[{c1==0,c2==0},{delv1,delv2}][[1]];
];

(* relation of vev and Lagrangian parameters (name it SMEFT$vev and keep vev
as variable) *)
SMEFT$VEV = vev0 + dvev // Simplify;
SMEFT$vev = vev;

(* invert mu-vev relation *)
SMEFT$muH = Simplify[ Sqrt[hlambda vev^2/2 - 3 Lam vev^4 MBScalarList["phi"]/4 -
                           Lam^2 vev^6 MBScalarList["phi8"]/2] // SMEFTExpOrder,
		           Assumptions -> vev > 0];
SMEFT$muH = CombineCommonPowers[SMEFT$muH, Lam, SMEFT$ExpansionOrder];

(* find kinetic terms *)
KinH = OptimizeIndex[ GetKineticTerms[ LHeff, FlavorExpand->False ] ];
(* change Hd^2H to -(dH)^2 *)
KinH = KinH /. F_ del[del[F_, Index[Lorentz, aa_]], Index[Lorentz, aa_]] ->
               - del[F, Index[Lorentz, aa]]^2;
KinH = KinH /. GP del[del[GPbar, Index[Lorentz, aa_]], Index[Lorentz, aa_]] ->
               - del[GP, Index[Lorentz, aa]] del[GPbar, Index[Lorentz, aa]];
KinH = KinH /. GPbar del[del[GP, Index[Lorentz, aa_]], Index[Lorentz, aa_]] ->
               - del[GP, Index[Lorentz, aa]] del[GPbar, Index[Lorentz, aa]];

(* find normalization of Higgs kinetic terms *)
tmp = 2 KinH /. G0->0 /. GP->0 /. del[H, aa_] -> 1;
SMEFT$Hnorm = (Hnorm /.  Normal[ Solve[tmp == 1, Hnorm, Reals][[2]] ]) // SMEFTExpOrder;

If [ SMEFT$RxiGaugeStatus,
  tmp = 2 KinH /. H->0 /. GP->0 /. del[G0, aa_] -> 1;
  SMEFT$G0norm = (G0norm /. Normal[ Solve[tmp == 1, G0norm, Reals][[2]] ]) // SMEFTExpOrder;

  tmp = KinH /. H->0 /. G0->0 /.
                del[GP, Index[Lorentz, aa_]] del[GPbar, Index[Lorentz, aa_]] -> 1;
  SMEFT$GPnorm = (GPnorm /. Normal[ Solve[tmp == 1, GPnorm, Reals][[2]] ]) // SMEFTExpOrder;
,
  SMEFT$G0norm = 1;
  SMEFT$GPnorm = 1;
];

SMEFT$Hnorm  = 1/(1/SMEFT$Hnorm  // SMEFTExpOrder);
SMEFT$G0norm = 1/(1/SMEFT$G0norm // SMEFTExpOrder);
SMEFT$GPnorm = 1/(1/SMEFT$GPnorm // SMEFTExpOrder);

(* Higgs Lagrangian in normalized fields H,G *)
LHeff = LHeff // Expand;
LHeff = LHeff // SMEFTExpOrder;

(* find physical Higgs mass^2 *)
VGHiggs = VGHiggs /. Hnorm -> SMEFT$Hnorm // Expand // SMEFTExpOrder;
VGHiggs = VGHiggs /. vev -> SMEFT$VEV // Expand // SMEFTExpOrder;
MGH2 = GetMassTerms[VGHiggs]  // Expand // Simplify;

SMEFT$MH2 = 2 Coefficient[MGH2, H^2] /. muH -> SMEFT$muH // Expand // SMEFTExpOrder;
SMEFT$MH  = SMEFTExpOrder[ Sqrt[SMEFT$MH2] ];
SMEFT$MH  = Simplify[SMEFT$MH, Assumptions -> vev>0];
SMEFT$MH  = CombineCommonPowers[SMEFT$MH, Lam, SMEFT$ExpansionOrder];

If [ SMEFT$RxiGaugeStatus,
  MG0 = 2 Coefficient[MGH2, G0^2] /. G0norm->SMEFT$G0norm // Expand // SMEFTExpOrder ;
  MGP = Coefficient[MGH2, GP GPbar] /. GPnorm->SMEFT$GPnorm // Expand // SMEFTExpOrder;
  If[ MG0 =!= 0 || MGP =!= 0, Print["Error: non-vanishing Goldstone masses in unitary gauge!"]; Abort[];];
];

]]
(* end of HiggsSectorDiagonalization *)
]




GaugeSectorDiagonalization = Function[{},
(* GAUGE SECTOR REDIAGONALIZATION *)
Module[{aa, bb, cc, dd, k, mu$1, mu$2, tmp, LGeff, KinG, KinG0, KinF,
KinZ, KinZF, KinW, MassG, MassG0, MassW, MassF, MGKin, MGMass, MGKin1,
MGMass1, CKF, CKZ, CKZF, CMF, CMZ, CMZF, U0, U1, g1, gw},

(* Z and photon first

   | W3 |          | Z |
   |    | = AZnorm |   |
   | B  |          | A |

*)

(* pure gauge Lagrangian, X Xtilde neglected -> full derivatives in bilinears *)
LGeff = SMEFT$LGeff /. H -> 0 /. G0 -> 0 /. GP->0 /. GPbar->0 // Expand;
LGeff = LGeff /. Eps[___]-> 0;

KinG = GetKineticTerms[ LGeff, FlavorExpand->False ];
MassG = GetMassTerms[ LGeff, FlavorExpand->False ];

(* kinetic and mass terms for B and W3 fields *)
KinG0 = KinG /. W[___]-> 0 /. Wbar[___]-> 0 /. G[___] -> 0;
MassG0 = MassG /. W[___] -> 0 /. Wbar[___]-> 0 /. G[___] -> 0;

KinG = KinG - KinG0;
MassG = MassG - MassG0;
(* optimize indices and antisymmetrize *)
KinG0 = OptimizeIndex[KinG0];
KinG0 = KinG0 /. del[cc_[aa_], bb_] -> 1/2 (del[cc[aa], bb] - del[cc[bb], aa]);
MassG0 = OptimizeIndex[MassG0];

KinF = KinG0 /. Z[___] -> 0 // Simplify;
KinZ = KinG0 /. A[___] -> 0 // Simplify;
KinZF = KinG0 - KinZ - KinF // Expand // Simplify;

(* coefficients of kinetic terms *)
CKF =  -4 (KinF /. del[aa__] - del[bb__] -> 1);
CKZ =  -4 (KinZ /. del[aa__] - del[bb__] -> 1);
CKZF = -4 (KinZF /. del[aa__] - del[bb__] -> 1);

(* kinetic mixing matrix, kinetic term is -1/4 (dW3,dB) MGKin (dW3,dB) *)
MGKin = {{CKZ,CKZF/2},{CKZF/2,CKF}};

(* coefficients of mass terms *)
CMF = MassG0 /. Z[__]->0 // Simplify;
CMZ = MassG0 /. A[__]->0 // Simplify;
CMZF = MassG0 - CMF - CMZ // Simplify;

(* kinetic mixing matrix, kinetic term is 1/2 (B,W3) MGMAss (B,W3) *)
MGMass = 2 {{CMZ,CMZF/2},{CMZF/2,CMF}} /. Z[__] -> 1 /. A[__] -> 1;

(* AZnorm should be SM rotation times higher 1/Lambda^2 correction *)

g1 = g1norm G1;
gw = gwnorm GW;

(* First rotate out 1/Lambda^2 terms *)
U0 = 1/Sqrt[gw^2 + g1^2] {{gw, g1},{-g1, gw}};

MGKin1 = MGKin;
MGMass1 = MGMass;
 For[aa=1,aa<3,aa++,
  For[bb=1,bb<3,bb++,
    MGKin1 = MGKin1 /. AZnorm[Index[SU2W, aa], Index[SU2W, bb]] -> U0[[aa, bb]] // Expand;
    MGMass1 = MGMass1 /. AZnorm[Index[SU2W, aa], Index[SU2W, bb]] -> U0[[aa, bb]] // Expand;
  ];
];

MGKin1 = MGKin1 // FullSimplify;
MGMass1 = MGMass1 // FullSimplify;

aa = MGKin1[[1,1]];
bb = MGKin1[[1,2]];
cc = MGKin1[[2,2]];

U1 = {{cc,0},{-bb, Sqrt[aa cc - bb^2]}};

SMEFT$AZnorm = U0.U1/Sqrt[cc (aa cc - bb^2)] // Simplify;

MGKin1 = MGKin;
MGMass1 = MGMass;
For[aa=1,aa<3,aa++,
  For[bb=1,bb<3,bb++,
    MGKin1 = MGKin1 /. AZnorm[Index[SU2W, aa], Index[SU2W, bb]] -> SMEFT$AZnorm[[aa, bb]] // Expand;
    MGMass1 = MGMass1 /. AZnorm[Index[SU2W, aa], Index[SU2W, bb]] -> SMEFT$AZnorm[[aa, bb]] // Expand;
  ];
];

MGKin1 = MGKin1 // Simplify;
MGMass1 = MGMass1 //Simplify;

If [ (MGKin1 =!= {{1,0},{0,1}}) || (MGMass1[[2,2]] =!=0), Print["Incorrect gauge sector diagonalization!"]; Abort[]; ];

(* Corrected Z mass and expanded AZ mixing*)
SMEFT$MZ2 = MGMass1[[1,1]] /. gwnorm -> SMEFT$gwnorm /. g1norm -> SMEFT$g1norm;
SMEFT$MZ2 = 4 SMEFTExpOrder[ SMEFT$MZ2 ]  / vev^2 / (G1^2+GW^2) // Expand;
SMEFT$MZ = SMEFTExpOrder[ Sqrt[SMEFT$MZ2] ] // Expand;
SMEFT$MZ2 = (G1^2+GW^2) vev^2/4 CombineCommonPowers[SMEFT$MZ2,Lam,SMEFT$ExpansionOrder];
SMEFT$MZ = Sqrt[G1^2+GW^2] vev/2 CombineCommonPowers[SMEFT$MZ,Lam,SMEFT$ExpansionOrder];

SMEFT$AZnormNoExp = SMEFT$AZnorm;
SMEFT$AZnorm = SMEFT$AZnorm /. gwnorm -> SMEFT$gwnorm /. g1norm -> SMEFT$g1norm;
SMEFT$AZnorm = Sqrt[G1^2 + GW^2] SMEFTExpOrder[SMEFT$AZnorm] // Expand;
SMEFT$AZnorm = CombineCommonPowers[SMEFT$AZnorm,Lam,SMEFT$ExpansionOrder]/Sqrt[G1^2 + GW^2];

(* W+/W- sector rediagonalization, this is the easy one... *)
(* kinetic and mass terms for W+/W- fields *)
KinW = KinG /. G[___] -> 0 // Simplify;
MassW = MassG /. G[___] -> 0 // Simplify;
KinG = KinG - KinW // Expand // Simplify;
MassG = MassG - MassW // Expand // Simplify;

KinW = OptimizeIndex[KinW];
KinW = KinW /. del[cc_[aa_], bb_] -> 1/2 (del[cc[aa], bb] - del[cc[bb], aa]) // Simplify;
MassW = OptimizeIndex[MassW] // Simplify;

KinG = OptimizeIndex[KinG];
KinG = KinG /. del[cc_[aa_,dd_], bb_] -> 1/2 (del[cc[aa,dd], bb] - del[cc[bb,dd], aa]) // Simplify;

tmp = -2 (KinW /. del[aa__] - del[bb__] -> 1) // Simplify;

(* normalization of W field *)
tmp = Normal[Solve[tmp == 1, Wnorm, Reals][[2]]];
tmp = SMEFTExpOrder[Wnorm /. tmp];
SMEFT$Wnorm = CombineCommonPowers[tmp,Lam,SMEFT$ExpansionOrder];

(* corrected W mass *)
SMEFT$MW2 = MassW /. W[__]->1 /. Wbar[__]->1 /. gwnorm->SMEFT$gwnorm /. Wnorm->SMEFT$Wnorm // SMEFTExpOrder;
SMEFT$MW2 = CombineCommonPowers[SMEFT$MW2,Lam,SMEFT$ExpansionOrder];
SMEFT$MW  = SMEFTExpOrder[Sqrt[SMEFT$MW2]] /. Sqrt[GW^2 vev^2]-> GW vev;
SMEFT$MW  = CombineCommonPowers[SMEFT$MW,Lam,SMEFT$ExpansionOrder];

(* Gluon sector rediagonalization, again easy one... *)
(* normalization of gluon field *)
tmp = -4 (KinG /. del[aa__] - del[vb__] -> 1) // Simplify;
tmp = Normal[Solve[tmp == 1, Gnorm, Reals][[2]]];
tmp = SMEFTExpOrder[Gnorm /. tmp];
SMEFT$Gnorm = CombineCommonPowers[tmp,Lam,SMEFT$ExpansionOrder];

]
(* end of GaugeSectorDiagonalization *)
];





FermionSectorDiagonalization = Function[{},
(* diagonalization of the fermion sector *)
Module[{aa, bb, cc, k, LFeff, KinF, KinF0, MassF, Flist, MLL, MUU, MDD, lhs, sol,
xl, xd, xu, f1, f2, g1, g2},

KinF0 = I IndexDelta[Index[Generation, Generation$1], Index[Generation, Generation$2]]
        ((dqbar[SP$1, Generation$1, Colour$1].del[dq[SP$2, Generation$2, Colour$1], mu$1] +
          lbar[SP$1, Generation$1].del[l[SP$2, Generation$2], mu$1] +
	  uqbar[SP$1, Generation$1, Colour$1].del[uq[SP$2, Generation$2, Colour$1], mu$1]) *
	   Ga[mu$1, SP$1, SP$2] +
	  vlbar[SP$1, Generation$1].del[vl[SP$2, Generation$2], mu$1] *
	   TensDot[Ga[mu$1], ProjM][SP$1, SP$2]) // OptimizeIndex;

LFeff = SMEFT$LGFmass /. Z[___] -> 0 /. A[___] -> 0 /. W[___] -> 0 // Expand;
LFeff = LFeff // SMEFTExpOrder;

(* apply unitarity of fermion rotations *)
LFeff = LFeff //. VVL[aa_,bb_] Conjugate[VVL[aa_,cc_]] -> IndexDelta[bb,cc];
LFeff = LFeff //. VLL[aa_,bb_] Conjugate[VLL[aa_,cc_]] -> IndexDelta[bb,cc];
LFeff = LFeff //. VLR[aa_,bb_] Conjugate[VLR[aa_,cc_]] -> IndexDelta[bb,cc];
LFeff = LFeff //. VUL[aa_,bb_] Conjugate[VUL[aa_,cc_]] -> IndexDelta[bb,cc];
LFeff = LFeff //. VUR[aa_,bb_] Conjugate[VUR[aa_,cc_]] -> IndexDelta[bb,cc];
LFeff = LFeff //. VDL[aa_,bb_] Conjugate[VDL[aa_,cc_]] -> IndexDelta[bb,cc];
LFeff = LFeff //. VDR[aa_,bb_] Conjugate[VDR[aa_,cc_]] -> IndexDelta[bb,cc];

(* kinetic terms already in canonical form, no corrections *)
KinF = OptimizeIndex[ GetKineticTerms[ LFeff, FlavorExpand->False ] ];
If[Expand[KinF - KinF0] !=0, Print["Non-canonical fermion kinetic terms, feature not implemented!"];Abort[]];
(* find new mass terms *)
MassF = OptimizeIndex[GetMassTerms[ LFeff, FlavorExpand->False ]];

(* use FeynmanRules to simplify indices *)
Flist = FeynmanRules[ H MassF, ScreenOutput->False ];
(* keep only the right part *)
Flist = Flist /. ProjM[___] -> 0;

For[aa=1,aa<4,aa++,
  If[ Flist[[aa,1,2,1]] === l,  MLL =  - I Flist[[aa,2]] // Expand ];
  If[ Flist[[aa,1,2,1]] === uq, MUU =  - I Flist[[aa,2]] // Expand ];
  If[ Flist[[aa,1,2,1]] === dq, MDD =  - I Flist[[aa,2]] // Expand ];
];

(* neglect Kronecker delta in color indices and chiral projectors *)
MDD = MDD /. ProjP[___] -> 1 /. IndexDelta[aa_,bb_] -> 1;
MUU = MUU /. ProjP[___] -> 1 /. IndexDelta[aa_,bb_] -> 1;
MLL = MLL /. ProjP[___] -> 1;

SMEFT$MDD = MDD;
SMEFT$MUU = MUU;
SMEFT$MLL = MLL;

(* find relation of Yukawas to VL,VR from the condition of diagonal
tree level fermion masses *)

(* another trick - VL, VR are unitary, so

VL[a,b] VL[a,c]^* = delta_ac
VL[b,a] VL[c,a]^* = delta_ac

Thus for invertions I take

1/VLL[Generation$1, f1] -> Conjugate[VLL[Generation$1,f1]]

etc: seems to work...

*)

f1 = Index[ Generation, Ext[1] ];
f2 = Index[ Generation, Ext[2] ];

lhs = MLL /. yl[aa_, bb_] -> xl /. VLL[aa_, bb_] -> 1/Conjugate[VLL[aa,bb]] /.
                                   VLR[aa_, bb_] -> 1/Conjugate[VLR[aa,bb]];

sol = xl /. Solve[ lhs == - fml[f1,f2], xl ][[1]] // SMEFTExpOrder;

SMEFT$YL[g1_,g2_] = sol /. Ext[1] -> Unique[YL] /. Ext[2] -> Unique[YL] /.
                   Index[Generation,Generation$1] -> g1 /.
                   Index[Generation,Generation$2] -> g2 // Simplify;

lhs = MDD /. yd[aa_, bb_] -> xd /. VDL[aa_, bb_] -> 1/Conjugate[VDL[aa,bb]] /.
                                   VDR[aa_, bb_] -> 1/Conjugate[VDR[aa,bb]];

sol = xd /. Solve[ lhs == - fmd[f1,f2], xd ][[1]] // SMEFTExpOrder;

SMEFT$YD[g1_,g2_] = sol /. Ext[1] -> Unique[YD] /. Ext[2] -> Unique[YD] /.
                   Index[Generation,Generation$1] -> g1 /.
                   Index[Generation,Generation$2] -> g2 // Simplify;


lhs = MUU /. yu[aa_, bb_] -> xu /. VUL[aa_, bb_] -> 1/Conjugate[VUL[aa,bb]] /.
                                   VUR[aa_, bb_] -> 1/Conjugate[VUR[aa,bb]];

sol = xu /. Solve[ lhs == - fmu[f1,f2], xu ][[1]] // SMEFTExpOrder;

SMEFT$YU[g1_,g2_] = sol /. Ext[1] -> Unique[YU] /. Ext[2] -> Unique[YU] /.
                   Index[Generation,Generation$1] -> g1 /.
                   Index[Generation,Generation$2] -> g2 // Simplify;

]
(* end of FermionSectorDiagonalization *)
];

