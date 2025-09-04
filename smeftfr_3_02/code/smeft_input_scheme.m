(* SmeftFR v3.0 package *)
(* Defines numerical values of input parameters and their relations to
   SM Lagrangian parameters

This file can be modified by users so they can change it to their
preferred input scheme. Currently, as a default choice, it uses
physical particle masses as the input parameters (other choices here
are technically very inconvenient as Wilson coefficients may appear in
denominators of propagators).  In addition it uses as input G_fermi
(calculated from the muon lifetime), alpha_s(MZ), Wolfenstein
parameters and PMNS mixing angles. Parameters in the electroweak and
Higgs sector of the SM Lagrangian are expressed in terms of input
observables consistently including 1/Lambda^4 terms. Strong coupling
constant is related to alpha_s(MZ) including only the factor giving
the canonical normalization of the gluon fields, other types of dim-6
and dim-8 corrections are not easily available in the literature.
Quark sector parameters (CKM matrix) are calculated including dim6
corrections as proposed in 1812.08163. Finally, PMNS mixing is
calculated using the SM relations, no EFT terms included as of now.
*)


(* *************************************** *)
(* ***** "GF input scheme"           ***** *)
(* *************************************** *)

(* add here all user-defined parameter lists and do NOT REMOVE this line! *)
SM$InputParameters := Join[SM$GaugeParameters, SM$FermionParameters, SM$AuxiliaryParameters];

(* all the user defined parameters should be real scalars, no indices! *)

SM$GaugeParameters = {
(* observables used as input parameters in gauge and Higgs sector *)
  alphas == {
    ParameterType    -> External,
    Value            -> 0.1176,
    InteractionOrder -> {QCD,2},
    TeX              -> Subscript[\[Alpha],s],
    Description      -> "average alpha_s at MZ scale"
   },
(* either Gf or AEM are used depending on the choice of input scheme, not both of them together! *)
  Gf == {
    ParameterType    -> External,
    Value            -> Sqrt[ 192 Pi^3/taumu/mumass^5 ], (* cannot be a formulae so updated later in the code *)
    InteractionOrder -> {QED,2},
    TeX              -> Subscript[G,F],
    Description      -> "G_Fermi"
	},
  AEM == {
    ParameterType    -> External,
    Value            -> 1/137.035999084, 
    InteractionOrder -> {QED,2},
    TeX              -> Subscript[\[Alpha],em],
    Description      -> "alpha_em at p^2=0"
	},
  Zmass == {
    ParameterType    -> External,
    Value            -> 91.1876,
    TeX              -> Subscript[M,Z],
    Description      -> "Pole Z boson mass"
  },
  Wmass == {
    ParameterType    -> External,
    Value            -> 80.379,
    TeX              -> Subscript[M,W],
    Description      -> "Pole W boson mass"
  },
  Hmass == {
    ParameterType    -> External,
    Value            -> 125.35,
    TeX              -> Subscript[M,H],
    Description      -> "Pole Higgs boson mass"
  }
};


SM$FermionParameters = {
(* observables used as input parameters in fermion sector*)
  tmass == {
    ParameterType    -> External,
    Value            -> 172.76,
    TeX              -> Subscript[m,t],
    Description      -> "Pole top quark mass"
  },
  cmass == {
    ParameterType    -> External,
    Value            -> 1.27,
    TeX              -> Subscript[m,c],
    Description      -> "c quark mass"
  },
  umass == {
    ParameterType    -> External,
    Value            -> 2.55 10^-3,
    TeX              -> Subscript[m,u],
    Description      -> "u quark mass"
  },
  bmass == {
    ParameterType    -> External,
    Value            -> 4.7,
    TeX              -> Subscript[m,b],
    Description      -> "b quark mass"
  },
  smass == {
    ParameterType    -> External,
    Value            -> 0.101,
    TeX              -> Subscript[m,s],
    Description      -> "s quark mass"
  },
  dmass == {
    ParameterType    -> External,
    Value            -> 5.04 10^-3,
    TeX              -> Subscript[m,d],
    Description      -> "d quark mass"
  },
  taumass == {
    ParameterType    -> External,
    Value            -> 1.77686,
    TeX              -> Subscript[m,\[Tau]],
    Description      -> "tau lepton mass"
  },
  mumass == {
    ParameterType    -> External,
    Value            -> 0.1056583755,
    TeX              -> Subscript[m,\[Mu]],
    Description      -> "muon mass"
  },
  emass == {
    ParameterType    -> External,
    Value            -> 0.00051099895000,
    TeX              -> Subscript[m,e],
    Description      -> "electron mass"
  },
(* Wolfenstein parameters from Table 3 of 1812.08163 *)
  lamT == {
    ParameterType    -> External,
    Value            -> 0.22537,
    TeX              -> Subscript[\[Lambda],T],
    Description      -> "Wolfenstein lambda parameter"
  },
  AT == {
    ParameterType    -> External,
    Value            -> 0.828,
    TeX              -> Subscript[A,T],
    Description      -> "Wolfenstein A parameter"
  },
  rhoT == {
    ParameterType    -> External,
    Value            -> 0.194,
    TeX              -> Subscript[\[Rho],T],
    Description      -> "Wolfenstein rho parameter"
  },
  etaT == {
    ParameterType    -> External,
    Value            -> 0.391,
    TeX              -> Subscript[\[Eta],T],
    Description      -> "Wolfenstein eta parameter"
  },
(* PMNS angles and phase *)
  theta12 == {
    ParameterType    -> External,
    Value            -> 33.44/180 Pi,
    TeX              -> Subscript[\[Theta],12],
    Description      -> "PMNS theta12 angle"
  },
  theta13 == {
    ParameterType    -> External,
    Value            -> 8.57/180 Pi,
    TeX              -> Subscript[\[Theta],13],
    Description      -> "PMNS theta13 angle"
  },
  theta23 == {
    ParameterType    -> External,
    Value            -> 49.2/180 Pi,
    TeX              -> Subscript[\[Theta],23],
    Description      -> "PMNS theta23 angle"
  },
  delta == {
    ParameterType    -> External,
    Value            -> 194.001/180 Pi,
    TeX              -> \[Delta],
    Description      -> "PMNS Dirac phase"
  }
};


SM$AuxiliaryParameters = {
(* various additional input observables *)
  taumu == {
    ParameterType    -> External,
(* 1 GeV^−1 = 6.582 10^−25 s, tau_mu in s: 2.1969811(22)×10−6 *)
    Value            -> 3.337862504 10^18,
    TeX              -> Subscript[\[Tau],\[Mu]],
    Description      -> "muon lifetime in GeV^-1"
  },
  mB == {
    ParameterType    -> External,
    Value            -> 5.27932,
    TeX              -> Superscript[m,"B+"],
    Description      -> "B^+ meson mass"
   },
  mBs == {
    ParameterType    -> External,
    Value            -> 5.36689,
    TeX              -> Subscript[m,Bs],
    Description      -> "B_s meson mass"
   },
  mBd == {
    ParameterType    -> External,
    Value            -> 5.27963,
    TeX              -> Subscript[m,Bd],
    Description      -> "B_d meson mass"
   },
  mK == {
    ParameterType    -> External,
    Value            -> 0.493677,
    TeX              -> Superscript[m,"K+"],
    Description      -> "K^+ meson mass"
   },
  mPi == {
    ParameterType    -> External,
    Value            -> 0.13957061,
    TeX              -> Superscript[m,"\[Pi]+"],
    Description      -> "pi^+ meson mass"
   },
  S1 == {
    ParameterType    -> External,
    Value            -> 2.3124,
    TeX              -> Subscript[S,1][Subscript[m,W]],
    Description      -> "S_1(m_W) as of Table 2 in 1812.08163"
   }
};

(* PART OF PARAMETER DEFINITION!
all input parameters must have numerical values assigned, not the
formulae, so update Gf with taumu and taumass here! Alternatively
just edit parameters above directly giving numerical value to Gf *)
Module[{numf,numa,x},
numf = Rule @@@ Partition[ Riffle[ SM$FermionParameters[[All, 1]],
	  Value /. SM$FermionParameters[[All, 2]]], 2];
numa = Rule @@@ Partition[ Riffle[ SM$AuxiliaryParameters[[All, 1]],
          Value /. SM$AuxiliaryParameters[[All, 2]]], 2];

x = Position[SM$GaugeParameters[[All, 1]], Gf][[1, 1]];
SM$GaugeParameters[[x]] = SM$GaugeParameters[[x]] /. numf /. numa;
];



SMEFTInputScheme[ OptionsPattern[ { InputScheme -> "GF"} ] ] :=
(* routine calling initialization of parameters in various sectors -
gauge, quark and lepton. Has to be called BEFORE starting FeynRules,
reading in model files and calculating bilinears! *)
Module[{},

(* note that corrections from WCs to CKM matrix are always calculated
but added to only if SMEFTInitializeModel was called with option
CKMInput-> "yes" or "force" *)
   
Switch[ ToUpperCase[ OptionValue[InputScheme] ],
  "GF",  SMEFTGaugeInputsGF[]; 
         SMEFTQuarkInputs[]; 
         SMEFTLeptonInputs[],
  "AEM", SMEFTGaugeInputsAEM[];
         SMEFTQuarkInputs[];
         SMEFTLeptonInputs[],
   _,    Print["Only GF and AEM input shemes currently predefined in SmeftFR v3."];
         Print["To use scheme " <> OptionValue[InputScheme] <> " users must provide their own routine " <>
	       "relating default SMEFT parameters with their preferred set of input parameters," <>
	       " see manual for more explanations. Currently undefined, aborting..."];
         Abort[];
];

]
(* end of SMEFTInputScheme *)

				



SMEFTGaugeInputsGF = Function[{},
(* Gf input scheme:

in this function user has to define relations between physical
obsevables used as input parameters and the set of gauge and Higgs
related variables used to calculate Feynman rules by SmeftFR.
Parameters which needs to be defined here are:

   G,G',alpha_s (effective gauge couplings, see manual of SmeftFR)
   v, lambda (Higgs vev and quartic coupling)
   MW,MZ,MH (physical W,Z and H masses)

Expressions are expanded in 1/Lambda powers to SMEFT$$ExpansionOrder
   (0-> SM, 1->Dim6, 2->Dim8+Dim6^2) *)
			
Block[{Gf, MH, MW, MZ, Hmass, Zmass, Wmass}, Module[{aa, bb, cc, dd,
l1, l3, ll, le, c6, c8, V, Zgb, Zgw, Zgw3, Zg0, Zh, eps, g, gp,
hlambda},

Print["\nUsing (G_F, MZ, MW, MH) input scheme for the EW sector... "];
			       
$Assumptions = Gf>0 && Gf\[Element]Reals && MZ>MW;

(* default SMEFT parameters evaluation *)
MH = Hmass;
MZ = Zmass;
MW = Wmass;
      
V = Sqrt[(Sqrt[4 (2 Gf^2 - Lam^2 c8) + Lam^2 c6^2] + Lam c6)/2/(2 Gf^2 - Lam^2 c8)] // SMEFTExpOrder;
		      
Zgb = Sqrt[1 - 2 Lam V^2 MBScalarList["phiB"] - Lam^2 V^4 MBScalarList["B2phi4n1"] ] // SMEFTExpOrder;
Zgw = Sqrt[1 - 2 Lam V^2 MBScalarList["phiW"] - Lam^2 V^4 MBScalarList["W2phi4n1"] ] // SMEFTExpOrder;
Zg0 = Sqrt[1 + Lam/2 V^2 MBScalarList["phiD"] + Lam^2/4 V^4 MBScalarList["phi6D2"] ] // SMEFTExpOrder;
Zgw3 = Zgw - Lam^2 V^4/2 MBScalarList["W2phi4n3"] // SMEFTExpOrder;
(* Zgw3 = Sqrt[Zgw^2 - Lam^2 V^4 MBScalarList["W2phi4n3"]] // SMEFTExpOrder;*)
			
g   = 2 MW/V // SMEFTExpOrder;
eps = Lam V^2/Zgw3/Zgb (MBScalarList["phiWB"] + Lam V^2/2 MBScalarList["WBphi4n1"]) // SMEFTExpOrder;
gp  = Sqrt[((2 MZ/V/Zg0)^2 - (g Zgw/Zgw3)^2) (1 - eps^2)] - eps g Zgw/Zgw3 // SMEFTExpOrder // Simplify;
gp  = Simplify[gp /. Sqrt[Gf aa_] -> Sqrt[Gf] Sqrt[aa]];

Zh  = Sqrt[1 + Lam/2 V^2 MBScalarList["phiD"] + Lam^2/4 V^4 MBScalarList["phi6D2"] -
               2 Lam V^2 MBScalarList["phiBox"] - Lam^2 V^4 MBScalarList["phi6Box"]] // SMEFTExpOrder;
hlambda = Zh^2 MH^2/V^2 + 3 Lam V^2 MBScalarList["phi"] + 3 Lam^2 V^4 MBScalarList["phi8"] // SMEFTExpOrder;
			
l1[aa_,bb_] = MBTensor2List["phil1",aa,bb];
l3[aa_,bb_] = MBTensor2List["phil3",aa,bb];
ll[aa_,bb_,cc_,dd_] = MBTensor4List["ll",aa,bb,cc,dd];
le[aa_,bb_,cc_,dd_] = MBTensor4List["le",aa,bb,cc,dd];

c6 = -2 (ll[2,1,1,2] - l3[1,1] - l3[2,2]) // ExpandIndices;

c8 = ll[2,1,1,2]^2 - 2 ll[2,1,1,2] l3[1,1] - 2 ll[2,1,1,2] l3[2,2] +
     l1[1,2] l1[2,1] - l3[1,+2] l3[2,1] + l1[2,1] l3[1,2] - l1[1,2]
     l3[2,1] + l3[1,1]^2 + l3[2,2]^2 + 4 l3[1,1] l3[2,2] +
     le[2,1,1,2]^2/4 // ExpandIndices;

(* store analytical expanded expressions *)			
UserInput$vev  = V // SMEFTExpOrder;
UserInput$GW   = g // SMEFTExpOrder;
UserInput$G1   = gp // SMEFTExpOrder;
UserInput$GS   = 2 Sqrt[Pi] Sqrt[alphas] // SMEFTExpOrder;
UserInput$hlambda = hlambda // SMEFTExpOrder;
UserInput$MZ   = Zmass // SMEFTExpOrder;
UserInput$MW   = Wmass // SMEFTExpOrder;
UserInput$MH   = Hmass // SMEFTExpOrder;

]]
(* end of SMEFTGaugeInputs *)
]





SMEFTGaugeInputsAEM = Function[{},
(* alpha_em input scheme:

in this function user has to define relations between physical
obsevables used as input parameters and the set of gauge and Higgs
related variables used to calculate Feynman rules by SmeftFR.
Parameters which needs to be defined here are:

   G,G',G_s (effective gauge couplings, see manual of SmeftFR)
   v, lambda (Higgs vev and quartic coupling)
   MW,MZ,MH (physical W,Z and H masses)

Expressions are expanded in 1/Lambda powers to SMEFT$ExpansionOrder
   (0-> SM, 1->Dim6, 2->Dim8+Dim6^2) *)

Block[{AEM, MH, MW, MZ, Hmass, Zmass, Wmass}, Module[{gSM, gD6, gD8,
gpSM, gpD6, gpD8, vSM, vD6, vD8, Zh, V, g, gp, hlambda},

Print["\nUsing (alpha_em, MZ, MW, MH) input scheme for the EW sector... "];

MH = Hmass;
MZ = Zmass;
MW = Wmass;
				
gSM   = 2 MZ/Sqrt[MZ^2 - MW^2] Sqrt[Pi AEM];
gD6   = MW^3/(4 Pi AEM MZ^2) (MW MBScalarList["phiD"] + 4 Sqrt[MZ^2 - MW^2] MBScalarList["phiWB"]);
gD8   = MW^5/(32 Pi^2 AEM^2 MZ^4) (MZ^2 - MW^2) (- MW^3/(MZ^2 - MW^2) MBScalarList["phiD"]^2 +
      8 (MZ^2 - 3 MW^2)/Sqrt[MZ^2 - MW^2] MBScalarList["phiD"] MBScalarList["phiWB"] +
      32 (Sqrt[MZ^2 - MW^2] (MBScalarList["phiB"] + MBScalarList["phiW"]) - MW MBScalarList["phiWB"]) MBScalarList["phiWB"] +
      4 (MW MBScalarList["phi6D2"] + 8 MW MBScalarList["W2phi4n3"] -
      4 MZ^2/MW MBScalarList["W2phi4n3"] + 4 Sqrt[MZ^2 - MW^2] MBScalarList["WBphi4n1"]));

gpSM   = 2 MZ/MW Sqrt[Pi AEM];
gpD6   = -  MW^2 (MZ^2 - MW^2)/(4 Pi AEM MZ^2) MBScalarList["phiD"];
gpD8   = MW^4 (MZ^2 - MW^2)/(32 Pi^2 AEM^2 MZ^4) ((MW^2 + 3 MZ^2) MBScalarList["phiD"]^2 -
         16 (MZ^2 - MW^2)  MBScalarList["phiWB"]^2 +
         16 MW Sqrt[MZ^2 - MW^2] MBScalarList["phiD"] MBScalarList["phiWB"] -
         4 (MZ^2 - MW^2) (MBScalarList["phi6D2"] + 4 MBScalarList["W2phi4n3"]));

vSM   = MW/MZ Sqrt[MZ^2 - MW^2]/Sqrt[Pi AEM] ;
vD6   = - gD6;
vD8   = MW^5 (MZ^2 - MW^2)/(32 Pi^2 AEM^2 MZ^4) (3 MW^3/(MZ^2 - MW^2) MBScalarList["phiD"]^2 +
      8 (5 MW^2 - MZ^2)/Sqrt[MZ^2 - MW^2] MBScalarList["phiD"] MBScalarList["phiWB"] -    
     32 (Sqrt[MZ^2 - MW^2] (MBScalarList["phiB"] + MBScalarList["phiW"]) - 2 MW MBScalarList["phiWB"]) MBScalarList["phiWB"] -
     4 (MBScalarList["phi6D2"] MW + 4 (2 MW^2 - MZ^2)/MW MBScalarList["W2phi4n3"] + 4 Sqrt[MZ^2 - MW^2] MBScalarList["WBphi4n1"]));      

g   = gSM  (1 + Lam gD6  + Lam^2 gD8 ) // SMEFTExpOrder;
gp  = gpSM (1 + Lam gpD6 + Lam^2 gpD8) // SMEFTExpOrder;
V   = vSM  (1 + Lam vD6  + Lam^2 vD8 ) // SMEFTExpOrder;

Zh  = Sqrt[1 + Lam/2 V^2 MBScalarList["phiD"] + Lam^2/4 V^4 MBScalarList["phi6D2"] -
               2 Lam V^2 MBScalarList["phiBox"] - Lam^2 V^4 MBScalarList["phi6Box"]] // SMEFTExpOrder;
hlambda = Zh^2 MH^2/V^2 + 3 Lam V^2 MBScalarList["phi"] + 3 Lam^2 V^4 MBScalarList["phi8"] // SMEFTExpOrder;

(* store analytical expanded expressions *)
UserInput$vev  = V // SMEFTExpOrder;
UserInput$GW   = g // SMEFTExpOrder;
UserInput$G1   = gp // SMEFTExpOrder;
UserInput$GS   = 2 Sqrt[Pi] Sqrt[alphas] // SMEFTExpOrder;
UserInput$hlambda = hlambda // SMEFTExpOrder;
UserInput$MZ   = Zmass // SMEFTExpOrder;
UserInput$MW   = Wmass // SMEFTExpOrder;
UserInput$MH   = Hmass // SMEFTExpOrder;
      
]]
(* end of SMEFTGaugeInputsAEM *)
]




(* CAUTION: In what follows flavored operators are assumed to be in
MASS BASIS, i.e. including redefinitions to include rotations in
flavor space required to make fermion mass matrices diagonal *)

SMEFTQuarkInputs = Function[{},
(* in this function user has to define relations between physical
observables used as input parameters and quark sector masses and mixing.
Expressions are currently expanded to dimension 6 only, i.e. 1/Lambda^2

Current version of this routine uses the approach of 1812.08163
*)

Block[{Gf, AEM, MH, MW, MZ}, Module[{v, dvv, VCKM, VCKMT, VTW, VW,
LVLLvedu, LSRRvedu, LTRRvedu, LVLRvedu, LSRLvedu, epsA, epsP,
LVLLddSM, LVLLddSMa, LVLLdd, LVRRdd, LV1LRdd, LV8LRdd, C1, C1SM,
C1tilde, C4, C5, A, B1, B4, B5, a1, b, c, d, e, m, j, q, dlam, dA,
drho, deta, Dlam, DA, Drho, Deta, numg, numf, numa, DelW, eps, dV,
DV},

MZ = UserInput$MZ;
MW = UserInput$MW;
MH = UserInput$MH;
				    
$Assumptions = lamT > 0 && AT > 0 && rhoT > 0  && etaT > 0;
$Assumptions = $Assumptions && lamT\[Element]Reals && AT\[Element]Reals &&
                               rhoT\[Element]Reals && etaT\[Element]Reals;
$Assumptions = $Assumptions && Gf > 0 && Lam > 0 && Gf\[Element]Reals && Lam\[Element]Reals;

(* numerical values of parameters as substitution lists *)
numg = Rule @@@ Partition[ Riffle[ SM$GaugeParameters[[All, 1]],
	  Value /. SM$GaugeParameters[[All, 2]]], 2];
numf = Rule @@@ Partition[ Riffle[ SM$FermionParameters[[All, 1]],
	  Value /. SM$FermionParameters[[All, 2]]], 2];
numa = Rule @@@ Partition[ Riffle[ SM$AuxiliaryParameters[[All, 1]],
	  Value /. SM$AuxiliaryParameters[[All, 2]]], 2];

(* update input parameters for quark masses *)
UserInput$MQU = umass;
UserInput$MQC = cmass;
UserInput$MQT = tmass;
UserInput$MQD = dmass;
UserInput$MQS = smass;
UserInput$MQB = bmass;

(* corrections to vev at 1/Lambda^2 *)
v = UserInput$vev /. Lam->0;
dvv = Lam Coefficient[UserInput$vev, Lam];

(* VCKMT is "tilde" matrix containing numerical inputs, actual CKM
matrix is given by VCKM=VCKMT-dV, where dV is NP shift (see
e.g. 5.9). We need analytical (VTW) and numerical (VCKMT) form *)
VTW = Table[0,{i,1,3},{j,1,3}];
VTW[[1,1]] = 1 - lamT^2/2 - lamT^4/8;
VTW[[1,2]] = lamT;
VTW[[1,3]] = AT lamT^3 (1 + lamT^2/2) (rhoT - I etaT);
VTW[[2,1]] = - lamT + AT^2 lamT^5 (1/2 - rhoT - I etaT);
VTW[[2,2]] = 1 - lamT^2/2 - lamT^4/8 (1 + 4 AT^2);
VTW[[2,3]] = AT lamT^2 ;
VTW[[3,1]] = AT lamT^3 (1 - rhoT - I etaT);
VTW[[3,2]] = - AT lamT^2 + AT lamT^4 (1/2 - rhoT - I etaT);
VTW[[3,3]] = 1 - AT^2 lamT^4/2;

VCKMT = VTW /. numf;

(* Two semileptonic decay rates as experimental input *)

(* LEFT coefficients in terms of SMEFT at EW scale, eq. A.1 *)
LVLLvedu[i_,i_,x_,k_] := 2 Sum[ Conjugate[ VCKMT[[k,j]] ] MBTensor4List["lq3",i,i,x,j], {j,1,3} ] -
                         2 Conjugate[ Sum[ VCKMT[[k,j]] MBTensor2List["phiq3",j,x], {j,1,3} ] ] -
                         2 Conjugate[ VCKMT[[k,x]] ] MBTensor2List["phil3",i,i];

LSRRvedu[i_,i_,x_,k_] := MBTensor4List["lequ1",i,i,x,k];
LTRRvedu[i_,i_,x_,k_] := MBTensor4List["lequ3",i,i,x,k];
LVLRvedu[i_,i_,x_,k_] := - Conjugate[MBTensor2List["phiud",k,x]];
LSRLvedu[i_,i_,x_,k_] := Sum[ Conjugate[VCKMT[[k,j]] ] MBTensor4List["ledq",i,i,x,j], {j,1,3} ];

(* RGE evolution of eps parameters up to \[Mu]=MW, eq. B.1 *)
epsA[2,1,m_] := - v^2/2/VCKMT[[1,m]] Conjugate[1.0094 LVLLvedu[2,2,m,1] - 1.0047 LVLRvedu[2,2,m,1] ];
epsP[2,1,m_] := - v^2/2/VCKMT[[1,m]] Conjugate[1.73 (LSRRvedu[2,2,m,1] - LSRLvedu[2,2,m,1]) -
                                               0.024 LTRRvedu[2,2,m,1] ];

epsA[3,1,3] := - v^2/2/VCKMT[[1,3]] Conjugate[1.0075 LVLLvedu[3,3,3,1] - 1.0038 LVLRvedu[3,3,3,1] ];
epsP[3,1,3] := - v^2/2/VCKMT[[1,3]] Conjugate[1.45 (LSRRvedu[3,3,3,1] - LSRLvedu[3,3,3,1]) -
		                              0.018 LTRRvedu[3,3,3,1] ];

(* Equations 4.3*)
DBt2 = 2 Lam (Re[epsA[3,1,3]] - mB^2/(umass + bmass)/taumass Re[epsP[3,1,3]]) + 4 dvv/v;
(* and 4.7*)
DKPi = Lam (2 Re[epsA[2,1,2] - epsA[2,1,1]] -
	    2/mumass (mK^2 Re[epsP[2,1,2]]/(umass + smass) -
		      mPi^2 Re[epsP[2,1,1]]/(umass + dmass)));

(* Two mass differences as the experimental input *)
(* LEFT coefficients in terms of SMEFT at EW scale, eq. A.2 *)
LVLLdd[m_,3,m_,3]  := Lam (MBTensor4List["qq1",m,3,m,3] + MBTensor4List["qq3",m,3,m,3]);
LVRRdd[m_,3,m_,3]  := Lam MBTensor4List["dd",m,3,m,3];
LV1LRdd[m_,3,m_,3] := Lam MBTensor4List["qd1",m,3,m,3];
LV8LRdd[m_,3,m_,3] := Lam MBTensor4List["qd8",m,3,m,3];

(*Wilson coefficients determining mass differences, eq. B.2*)
C1SM[q_]    :=  - 0.858 MW^2/32/Pi^2/v^4 (VCKMT[[3,q]] Conjugate[VCKMT[[3,3]]])^2 S1;
C1[q_]      :=  0.858 LVLLdd[q,3,q,3]/C1SM[q];
C1tilde[q_] :=  0.858 LVRRdd[q,3,q,3]/C1SM[q];
C4[q_]      := (- 0.755 LV1LRdd[q,3,q,3] - 1.94 LV8LRdd[q,3,q,3])/C1SM[q];
C5[q_]      := (- 1.856 LV1LRdd[q,3,q,3] + 0.237 LV8LRdd[q,3,q,3])/C1SM[q];

(* Remaining constants, numerical values of bag parameters in MSbar scheme taken from ref [49] Table 2 *)
A  = {1,-5/8,1/8,3/4,1/4};
B1 = {0.85,0.86};
B4 = {0.95,0.93};
B5 = {1.47,1.57};

(* Finally, mass differences, eq. 4.11 *)
DMqd = 4 dvv/v + Re[C1[1] + C1tilde[1] +
		  (mBd/(bmass+dmass))^2 (A[[4]] B4[[1]]/B1[[1]] C4[1] +
					 A[[5]] B5[[1]]/B1[[1]] C5[1])];

DMqs = 4 dvv/v + Re[C1[2] + C1tilde[2] +
		  (mBs/(bmass+smass))^2 (A[[4]] B4[[2]]/B1[[2]] C4[2] +
					 A[[5]] B5[[2]]/B1[[2]] C5[2])];

(* All parameters that will define shift of Wolfenstein parameters *)
DelW = {DKPi,DBt2,DMqd,DMqs};
DelW = DelW /. numg /. numf /. numa // Expand;

(* Full matrix of eq. 4.18; here change in definition of e comparing
   to 4.19, to have correct numerical factors in matrix, in e a->AT *)

a1 = (1 - 2 rhoT)/2;
b  = (etaT^2 + (1 - rhoT)^2)/2;
c  = (etaT^2 + rhoT^2)/2;
d  = etaT^2 - rhoT^2 + rhoT;
e  = lamT^2 (1 - AT lamT^2);

m = {
    {lamT/2 - lamT^3/2, 0, 0, 0},
    {- AT + AT lamT^2 + c AT lamT^4, - c e AT, b e AT, AT/2 - a1 e AT},
    {a1 - b lamT^2 + c (5 - 4 rhoT)/2 lamT^4, c (1 - 2 e a1), - b (1 - 2 a1 e), a1 (1 - 2 a1 e)},
    {- d/2/etaT + b rhoT lamT^2/etaT - c (2 d + 3 (rhoT - 1))/2/etaT lamT^4,
       c/etaT (1 - rhoT + d e), b/etaT (rhoT - d e), - d/2/etaT (1 - 2 a1 e)}
    } /. numf;

(* Finally, eq. 4.17: shifts in Wolfenstein parameters in terms of experimental input and WC's *)
{Dlam, DA, Drho, Deta} = m.DelW // Expand;

(* Now, shifts of CKM matrix elements in terms of inputs *)

(* CKM matrix in terms of Wolfenstein parameters plus linear corrections in NP, VW = VTW - dVTW *)
VW = VTW /. lamT -> lamT - eps dlam /. AT -> AT - eps dA /.
            etaT -> etaT - eps deta /. rhoT -> rhoT - eps drho // Expand;
VW = Normal[VW + O[eps]^2] /. eps->1;
dV = VTW - VW // Expand;

(* Final substitutions to get shifts in CKM matrix elements, this is
THE final result, shifts of CKM matrix elements in terms of inputs -
Wolfenstein parameters masses etc and Wilsons *)

DV = dV /. dlam -> Dlam /. dA -> DA /. drho -> Drho /. deta -> Deta;

(* Kill terms of higher order in lamT, max order 6 now*)
DV = (Expand[DV /. lamT -> eps lamT] // Function[y,Normal[y + O[eps]^6]]) /. eps -> 1;
(* semi numerical approximation, with numerical values of WC coefficients *)
DV = DV /. numf /. numa // Expand;
(* kill small terms: 1000/Lambda^2->0 *)
(* DV = Chop[DV,1000]; *)

(* SMEFT CKM matrix elements it terms of input coefficients and WCs *)

UserInput$CKM = VCKMT;
UserInput$dCKM = - DV // Expand; 

If [ SMEFT$CKMInput != "no" && CKMCorrectionsSize[], 
  Print[ Style[ "Corrections to CKM elements from higher dimension operators" <>
	" most likely violate constraints from flavor processes measurements! " <>
	"Consider modifying values of flavor off-diagonal Wilson coefficients?", Red] ];
  If [ SMEFT$CKMInput == "force",
    Print[ Style[ "Option SMEFTInitializeModel[ ForceCKMInput -> True ] have been used" <>
	  " - proceeding with calculations in spite of unrealistic CKM values...", Red] ];
  ,
    Print[ Style[ "Aborting - modify flavor-violating WCs " <>
        "or to proceed in spite of unrealistic CKM values " <>
	"initialize calculations using option " <>
	"SMEFTInitializeModel[ ForceCKMInput -> True  ].", Red] ];
     Abort[];
     ];
   ];

UserInput$CKM = VCKMT - If[ SMEFT$CKMInput != "no", DV, 0 ] /. SMEFT$WCValues;  

]]
(* end of SMEFTQuarkInputsGF *)
]



CKMCorrectionsSize = Function[{}, 
Module[{dCKM, dCKM0, wc0, a},

(* list with all WC->0 and Lam->1 *)       
wc0 = (# -> 0 & /@ Complement[Keys[SMEFT$WCValues], {Lam}]);
wc0 = Join[wc0,{Lam->1}];

(* find numerical corrections to CKM from each WC separately,
neglecting vanishing contributions *)       
dCKM = ({ToString[Keys[#]], UserInput$dCKM} /. # & /@ SMEFT$WCValues) /. wc0;
dCKM = dCKM /. Complex[a_, 0.] -> a;
dCKM0 = {#[[1]], {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}} & /@ dCKM;
dCKM = Complement[dCKM, dCKM0];
       
(* find norm of deviations from CKM neglecting values smaller than threshold *)
dCKM = ({#[[1]], Chop[Norm[#[[2]]/UserInput$CKM], SMEFT$CKMTreshold]} & /@ dCKM);
dCKM0 = ({#[[1]], 0} & /@ dCKM);
dCKM = Complement[dCKM, dCKM0];
       
(* print warnings *)       
Print[Style["WARNING:",Red,Bold], " Wilson coefficient " <>
      #[[1]] <> "=", ToExpression[#[[1]]] /. SMEFT$WCValues,
      " leads to total ",
      Round[100 #[[2]]],
      "% correction to CKM matrix elements!"] & /@ dCKM;  
       
(* return True if large corrections found *)
! TrueQ[ dCKM == {} ] 
       
]
(* end of CKMCorrectionsSize *)			      
]



SMEFTLeptonInputs = Function[{},
(* in this function user has to define relations between physical
observables used as input parameters and lepton sector masses and
mixing.  Currently expressions have no EFT dependence - no corrections
to PMNS calculated as of now.  *)

Module[{vsq, numf, s12 ,c12, s13, c13, s23, c23, eph},

(* numerical values of parameters as substitution lists *)
numf = Rule @@@ Partition[ Riffle[ SM$FermionParameters[[All, 1]],
	  Value /. SM$FermionParameters[[All, 2]]], 2];

vsq = Sqrt[Lam] (UserInput$vev /. Lam->0)^2;

(* update numerical values of lepton and neutrino masses *)
UserInput$MLE = emass;
UserInput$MLM = mumass;
UserInput$MLT = taumass;
UserInput$MVE = vsq Abs[MBTensor2List["vv",1,1]];
UserInput$MVM = vsq Abs[MBTensor2List["vv",2,2]];
UserInput$MVT = vsq Abs[MBTensor2List["vv",3,3]];

c12 = Cos[theta12];
s12 = Sin[theta12];
c13 = Cos[theta13];
s13 = Sin[theta13];
c23 = Cos[theta23];
s23 = Sin[theta23];
eph = Exp[I delta];

(* PMNS given as numerical values *)
UserInput$PMNS = {{c12 c13, s12 c13, s13/eph},
        {-s12 c23 - c12 s23 s13 eph, c12 c23 - s12 s23 s13 eph, s23 c13},
        {s12 s23 - c12 c23 s13 eph, -c12 s23 - s12 c23 s13 eph, c23 c13}} /. numf;

UserInput$delPMNS = {{0,0,0},
		     {0,0,0},
		     {0,0,0}};

]
(* end of SMEFTLeptonInputs *)
];


