(* SmeftFR v3.0 package *)
(* Calculation of Feynman rules for 4-fermion and B+L violating vertices *)


Fierz4f = Function[{expr},
(* 4-fermion Fierzing routines *)
Module[{a,b,tmp,res,find,find1,find2,find3,find4,fierz,mu$1},

tmp = expr /. TensDot[___][Index[Spin, Ext[1]], Index[Spin, Ext[2]]] -> 0;
res = expr - tmp // Expand;

find1 = Unique[fierz];
find2 = Unique[fierz];
find3 = Unique[fierz];
find4 = Unique[fierz];

tmp = tmp /. TensDot[Ga[Index[Lorentz,a_]],b_][Index[Spin,Ext[1]],Index[Spin,Ext[4]]] *
             TensDot[Ga[Index[Lorentz,a_]],b_][Index[Spin,Ext[3]],Index[Spin,Ext[2]]] ->
             TensDot[Ga[Index[Lorentz,a]],b][Index[Spin,find1],Index[Spin,find2]] *
             TensDot[Ga[Index[Lorentz,a]],b][Index[Spin,find3],Index[Spin,find4]];

tmp = tmp /. TensDot[Ga[Index[Lorentz,a_]],ProjM][Index[Spin,Ext[1]],Index[Spin,Ext[4]]] *
             TensDot[Ga[Index[Lorentz,a_]],ProjP][Index[Spin,Ext[3]],Index[Spin,Ext[2]]] ->
             ( - 2 ProjP[Index[Spin,find1],Index[Spin,find2]] ProjM[Index[Spin,find3],Index[Spin,find4]] );

tmp = tmp /. TensDot[Ga[Index[Lorentz,a_]],ProjP][Index[Spin,Ext[1]],Index[Spin,Ext[4]]] *
             TensDot[Ga[Index[Lorentz,a_]],ProjM][Index[Spin,Ext[3]],Index[Spin,Ext[2]]] ->
             ( - 2 ProjM[Index[Spin,find1],Index[Spin,find2]] ProjP[Index[Spin,find3],Index[Spin,find4]] );

tmp = tmp /. ProjM[Index[Spin,Ext[1]],Index[Spin,Ext[4]]] ProjP[Index[Spin,Ext[3]],Index[Spin,Ext[2]]] ->
             ( - 1/2 TensDot[Ga[Index[Lorentz,mu$1]],ProjP][Index[Spin,find1],Index[Spin,find2]] *
                   TensDot[Ga[Index[Lorentz,mu$1]],ProjM][Index[Spin,find3],Index[Spin,find4]] );

tmp = tmp /. ProjP[Index[Spin,Ext[1]],Index[Spin,Ext[4]]] ProjM[Index[Spin,Ext[3]],Index[Spin,Ext[2]]] ->
             ( - 1/2 TensDot[Ga[Index[Lorentz,\[Mu]]],ProjM][Index[Spin,find1],Index[Spin,find2]] *
                   TensDot[Ga[Index[Lorentz,\[Mu]]],ProjP][Index[Spin,find3],Index[Spin,find4]] );

tmp = tmp /. ProjP[Index[Spin,Ext[1]],Index[Spin,Ext[4]]] ProjP[Index[Spin,Ext[3]],Index[Spin,Ext[2]]] ->
             ( - 1/2 ProjP[Index[Spin,find1],Index[Spin,find2]] ProjP[Index[Spin,find3],Index[Spin,find4]] -
               1/8 SIG[\[Mu],\[Nu],ProjP,find1,find2] SIG[\[Mu],\[Nu],ProjP,find3,find4] ) ;

tmp = tmp /. ProjM[Index[Spin,Ext[1]],Index[Spin,Ext[4]]] ProjM[Index[Spin,Ext[3]],Index[Spin,Ext[2]]] ->
             ( - 1/2 ProjM[Index[Spin,find1],Index[Spin,find2]] ProjM[Index[Spin,find3],Index[Spin,find4]] -
               1/8 SIG[\[Mu],\[Nu],ProjM,find1,find2] SIG[\[Mu],\[Nu],ProjM,find3,find4] ) ;

find = Unique[fierz];

tmp = tmp /. Ext[4] -> find;
tmp = tmp /. Ext[2] -> Ext[4];
tmp = tmp /. find   -> Ext[2];

tmp = tmp /. find1 -> Ext[1] /.
             find2 -> Ext[2] /.
             find3 -> Ext[3] /.
             find4 -> Ext[4];

tmp + res

]
(* end of Fierz4f *)
]




SymmetrizeNeutrino4Current = Function[{x},
(* symmetrize vertices with 4 fermions, 2 of them being Majorana neutrinos *)
Module[{a,swap,i,frule,symlist},

If [ Length[x]>0,
    
  symlist = {};

  For[i=1, i < Length[x]+1, i++,
    If[ SubsetQ[ Flatten[ Take[ x[[i]][[1]], All, 1] ], {vlbar,vl} ] &&
        (Flatten[ Take[ x[[i]][[1]], All, 1] ] =!= {vlbar,vl,vlbar,vl}),
      frule = x[[i]];
      Clear[swap];
      frule = frule /. Ext[3] -> swap;
      frule = frule /. Ext[4] -> Ext[3];
      frule = frule /. swap -> Ext[4];
      frule = frule /. TensDot[Ga[a_], ProjM][Index[Spin,Ext[4]],Index[Spin,Ext[3]]] ->
                     - TensDot[Ga[a],  ProjP][Index[Spin,Ext[3]],Index[Spin,Ext[4]]];
      frule = frule /. TensDot[Ga[a_], ProjP][Index[Spin,Ext[4]],Index[Spin,Ext[3]]] ->
                     - TensDot[Ga[a],  ProjM][Index[Spin,Ext[3]],Index[Spin,Ext[4]]];
      AppendTo[symlist,frule];
    ];
  ];

  MergeVertices[x,symlist]
,
  x
]
     
]
(* end of SymmetrizeNeutrino4Current *)
]




FourFermionInteractions = Function[{},
(* 4-fermion couplings *)
Module[{a, b, c, d, s1, s2, mu$2, mu$1, i, find, fierz, tmp, ia1, ia2, oplist, WC, k},

Print[Style["Calculating 4-fermion vertices...",Bold]];

oplist = Intersection[ SMEFT$OperatorList, {"ll","ee","le"}];
FourLeptonLagrangian = 0;
For [i=1, i<Length[oplist]+1, i++,
   FourLeptonLagrangian = FourLeptonLagrangian + Lam ToExpression[ "LQ"<>oplist[[i]] ];
];
FourLeptonLagrangian = ReplaceFlavorRotations[FourLeptonLagrangian] // Expand;
FourLeptonLagrangian = FourLeptonLagrangian /. Lam^k_->If[k>1,0,Lam^k];
FourLeptonVertices = FeynmanRules[FourLeptonLagrangian];
If [ SMEFT$MajoranaNeutrino,
   FourLeptonVertices = SymmetrizeNeutrino4Current[FourLeptonVertices];
   i = Position[ FourLeptonVertices[[All, 1]], {{vlbar, 1}, {vl, 2}, {vlbar, 3}, {vl,4}} ];
   If[ Length[i]>0, FourLeptonVertices[[i[[1,1]], 2]] = Neutrino4Vertex[] ];
 ];

oplist = Intersection[ SMEFT$OperatorList, {"lq1","lq3","eu","ed","lu","ld","qe","ledq","lequ1","lequ3"} ];
TwoQuarkTwoLeptonLagrangian = 0;
For [i=1, i<Length[oplist]+1, i++,
   TwoQuarkTwoLeptonLagrangian = TwoQuarkTwoLeptonLagrangian + Lam ToExpression[ "LQ"<>oplist[[i]] ];
];
TwoQuarkTwoLeptonLagrangian = SU2Simplify[TwoQuarkTwoLeptonLagrangian];
TwoQuarkTwoLeptonLagrangian = ReplaceFlavorRotations[TwoQuarkTwoLeptonLagrangian] // Expand;
TwoQuarkTwoLeptonLagrangian = TwoQuarkTwoLeptonLagrangian /. Lam^k_->If[k>1,0,Lam^k];
TwoQuarkTwoLeptonVertices = FeynmanRules[TwoQuarkTwoLeptonLagrangian];
If[ SMEFT$MajoranaNeutrino, TwoQuarkTwoLeptonVertices = SymmetrizeNeutrino4Current[TwoQuarkTwoLeptonVertices] ];

oplist = Intersection[ SMEFT$OperatorList, {"qq1","qq3","uu","dd","ud1","ud8","qu1","qu8","qd1","qd8","quqd1","quqd8"} ];
FourQuarkLagrangian = 0;
For [i=1, i<Length[oplist]+1, i++,
   FourQuarkLagrangian = FourQuarkLagrangian + Lam ToExpression[ "LQ"<>oplist[[i]] ];
];
FourQuarkLagrangian = SU2Simplify[FourQuarkLagrangian];
FourQuarkLagrangian = ReplaceFlavorRotations[FourQuarkLagrangian] // Expand;
FourQuarkLagrangian = FourQuarkLagrangian /. Lam^k_->If[k>1,0,Lam^k];
FourQuarkVertices = FeynmanRules[FourQuarkLagrangian];

(* store vertices in new variables for some simplifications *)
FourLeptonVerticesOrig = FourLeptonVertices;
FourQuarkVerticesOrig = FourQuarkVertices;
TwoQuarkTwoLeptonVerticesOrig = TwoQuarkTwoLeptonVertices;

TwoQuarkTwoLeptonVertices = TwoQuarkTwoLeptonVertices //. TensDot[Ga[a_],Ga[b_],c_][s1_,s2_] -> (ME[a,b] c[s1,s2] - I SIG[a,b,c,s1,s2]);
TwoQuarkTwoLeptonVertices = TwoQuarkTwoLeptonVertices /. SIG[Index[Lorentz,mu$2],Index[Lorentz,mu$1],c_,s1_,s2_] ->
             (- SIG[Index[Lorentz,mu$1],Index[Lorentz,mu$2],c,s1,s2]) // Expand;
TwoQuarkTwoLeptonVertices = TwoQuarkTwoLeptonVertices /. ME[Index[Lorentz,a_],Index[Lorentz,b_]] SIG[Index[Lorentz,a_],Index[Lorentz,b_],c_,s1_,s2_] -> 0;
TwoQuarkTwoLeptonVertices = TwoQuarkTwoLeptonVertices /. IndexDelta[Index[Colour, Ext[a_]], Index[Colour, Ext[b_]]] -> 1;


(* fierzing and index rotation, whenewer possible. Currently commented
out, can affect evanescent operators! *)

(*
For [ i=1, i < Length[FourLeptonVertices] + 1, i++,
    If [ FourLeptonVertices[[i,1]] === {{lbar,1},{l,2},{lbar,3},{l,4}} ||
	FourLeptonVertices[[i,1]] === {{vlbar,1},{vl,2},{vlbar,3},{vl,4}},
	FourLeptonVertices[[i,2]] = Fierz4f[ FourLeptonVertices[[i,2]] ] ];
    FourLeptonVertices[[i,2]] = Fierz4f[ FourLeptonVertices[[i,2]] ];
];

For [ i=1, i < Length[FourQuarkVertices] + 1, i++,
    If [ FourQuarkVertices[[i,1]] === {{uqbar,1},{uq,2},{uqbar,3},{uq,4}} ||
	FourQuarkVertices[[i,1]] === {{dqbar,1},{dq,2},{dqbar,3},{dq,4}},
	FourQuarkVertices[[i,2]] = Fierz4f[ FourQuarkVertices[[i,2]] ] ];
];
*)

(* "manual" symmetrization of indices in ll, ee operators *)

If[ MemberQ[ SMEFT$OperatorList, "ll" ],

WC = ToExpression[SMEFT$MB <> "ll"];

FourLeptonVertices = FourLeptonVertices /.
    WC @@ {Index[Generation, Generation$1], Index[Generation, Generation$2],
           Index[Generation, Ext[1]], Index[Generation, Ext[2]]} ->
    WC @@ {Index[Generation, Ext[1]], Index[Generation, Ext[2]],
           Index[Generation, Generation$1], Index[Generation, Generation$2]};

FourLeptonVertices = FourLeptonVertices /.
    WC @@ {Index[Generation, Ext[3]], Index[Generation, Ext[2]],
           Index[Generation, Ext[1]], Index[Generation, Ext[4]]} ->
    WC @@ {Index[Generation, Ext[1]], Index[Generation, Ext[4]],
           Index[Generation, Ext[3]], Index[Generation, Ext[2]]};

FourLeptonVertices = FourLeptonVertices /.
    WC @@ {Index[Generation, Ext[3]], Index[Generation, Ext[4]],
           Index[Generation, Ext[1]], Index[Generation, Ext[2]]} ->
    WC @@ {Index[Generation, Ext[1]], Index[Generation, Ext[2]],
           Index[Generation, Ext[3]], Index[Generation, Ext[4]]};
];

If[ MemberQ[ SMEFT$OperatorList, "ee" ],
 WC = ToExpression[SMEFT$MB <> "ee"];

FourLeptonVertices = FourLeptonVertices /. WC @@ {Index[Generation, Ext[3]], Index[Generation,
 Ext[2]], Index[Generation, Ext[1]], Index[Generation, Ext[4]]} ->
 WC @@ {Index[Generation, Ext[1]], Index[Generation, Ext[4]],
 Index[Generation, Ext[3]], Index[Generation, Ext[2]]};

FourLeptonVertices = FourLeptonVertices /. WC @@ {Index[Generation, Ext[3]], Index[Generation,
 Ext[4]], Index[Generation, Ext[1]], Index[Generation, Ext[2]]} ->
 WC @@ {Index[Generation, Ext[1]], Index[Generation, Ext[2]],
 Index[Generation, Ext[3]], Index[Generation, Ext[4]]};

];


(* "manual" symmetrization of indices in qq3, qq1, dd, uu operators *)

If[ MemberQ[ SMEFT$OperatorList, "qq1" ],
 WC = ToExpression[SMEFT$MB <> "qq1"];

FourQuarkVertices = FourQuarkVertices /. WC @@ {Index[Generation, Generation$2], Index[Generation,
   Generation$1], Index[Generation, Ext[1]], Index[Generation, Ext[2]]} ->
   WC @@ {Index[Generation, Ext[1]], Index[Generation, Ext[2]],
   Index[Generation, Generation$2], Index[Generation, Generation$1]};

FourQuarkVertices = FourQuarkVertices /. WC @@ {Index[Generation, Generation$3], Index[Generation,
   Generation$2], Index[Generation, Generation$4], Index[Generation,
   Generation$1] } -> WC @@ {Index[Generation, Generation$4 ],
   Index[Generation, Generation$1], Index[Generation, Generation$3],
   Index[Generation, Generation$2] };

FourQuarkVertices = FourQuarkVertices /. WC @@ {Index[Generation, Generation$4], Index[Generation,
   Generation$2], Index[Generation, Generation$3], Index[Generation,
   Generation$1] } -> WC @@ {Index[Generation, Generation$3 ],
   Index[Generation, Generation$1], Index[Generation, Generation$4],
   Index[Generation, Generation$2] };

FourQuarkVertices = FourQuarkVertices /. WC @@ {Index[Generation, Ext[3]], Index[Generation,
   Ext[2]], Index[Generation, Ext[1]], Index[Generation, Ext[4]]} ->
   WC @@ {Index[Generation, Ext[1]], Index[Generation, Ext[4]],
   Index[Generation, Ext[3]], Index[Generation, Ext[2]]};

FourQuarkVertices = FourQuarkVertices /. WC @@ {Index[Generation, Ext[3]], Index[Generation,
   Ext[4]], Index[Generation, Ext[1]], Index[Generation, Ext[2]]} ->
   WC @@ {Index[Generation, Ext[1]], Index[Generation, Ext[2]],
   Index[Generation, Ext[3]], Index[Generation, Ext[4]]};

];


If[ MemberQ[ SMEFT$OperatorList, "qq3" ],
  WC = ToExpression[SMEFT$MB <> "qq3"];

FourQuarkVertices = FourQuarkVertices /. WC @@ {Index[Generation, Generation$2], Index[Generation,
   Ext[2]], Index[Generation, Ext[1]], Index[Generation,
   Generation$1]} -> WC @@ {Index[Generation, Ext[1]], Index[Generation,
   Generation$1 ], Index[Generation, Generation$2], Index[Generation,
   Ext[2]]};

FourQuarkVertices = FourQuarkVertices /. WC @@ {Index[Generation, Generation$2], Index[Generation,
   Generation$1], Index[Generation, Ext[1]], Index[Generation, Ext[2]]}->
   WC @@ {Index[Generation, Ext[1]], Index[Generation, Ext[2]],
   Index[Generation, Generation$2], Index[Generation, Generation$1]};

FourQuarkVertices = FourQuarkVertices /. WC @@ {Index[Generation, Generation$3], Index[Generation,
   Generation$2], Index[Generation, Generation$4], Index[Generation,
   Generation$1]} -> WC @@ {Index[Generation, Generation$4 ],
   Index[Generation, Generation$1], Index[Generation, Generation$3],
   Index[Generation, Generation$2]};

FourQuarkVertices = FourQuarkVertices /. WC @@ {Index[Generation, Generation$4], Index[Generation,
   Generation$2], Index[Generation, Generation$3], Index[Generation,
   Generation$1] } -> WC @@ {Index[Generation, Generation$3 ],
   Index[Generation, Generation$1], Index[Generation, Generation$4],
   Index[Generation, Generation$2]};

FourQuarkVertices = FourQuarkVertices /. WC @@ {Index[Generation, Ext[3]], Index[Generation,
   Ext[2]], Index[Generation, Ext[1]], Index[Generation, Ext[4]]} ->
   WC @@ {Index[Generation, Ext[1]], Index[Generation, Ext[4]],
   Index[Generation, Ext[3]], Index[Generation, Ext[2]]};

FourQuarkVertices = FourQuarkVertices /. WC @@ {Index[Generation, Ext[3]], Index[Generation,
   Ext[4]], Index[Generation, Ext[1]], Index[Generation, Ext[2]]} ->
   WC @@ {Index[Generation, Ext[1]], Index[Generation, Ext[2]],
   Index[Generation, Ext[3]], Index[Generation, Ext[4]]};

];

If[ MemberQ[ SMEFT$OperatorList, "uu" ],
  WC = ToExpression[SMEFT$MB <> "uu"];

  FourQuarkVertices = FourQuarkVertices /. WC @@ {Index[Generation, Ext[3]], Index[Generation,
  Ext[2]], Index[Generation, Ext[1]], Index[Generation, Ext[4]]} ->
  WC @@ {Index[Generation, Ext[1]], Index[Generation, Ext[4]],
  Index[Generation, Ext[3]], Index[Generation, Ext[2]]};

  FourQuarkVertices = FourQuarkVertices /. WC @@ {Index[Generation, Ext[3]], Index[Generation,
  Ext[4]], Index[Generation, Ext[1]], Index[Generation, Ext[2]]} ->
  WC @@ {Index[Generation, Ext[1]], Index[Generation, Ext[2]],
  Index[Generation, Ext[3]], Index[Generation, Ext[4]]};

];

If[ MemberQ[ SMEFT$OperatorList, "dd" ],
  WC = ToExpression[SMEFT$MB <> "dd"];

  FourQuarkVertices = FourQuarkVertices /. WC @@ {Index[Generation, Ext[3]], Index[Generation,
  Ext[2]], Index[Generation, Ext[1]], Index[Generation, Ext[4]]} ->
  WC @@ {Index[Generation, Ext[1]], Index[Generation, Ext[4]],
  Index[Generation, Ext[3]], Index[Generation, Ext[2]]};

  FourQuarkVertices = FourQuarkVertices /. WC @@ {Index[Generation, Ext[3]], Index[Generation,
  Ext[4]], Index[Generation, Ext[1]], Index[Generation, Ext[2]]} ->
  WC @@ {Index[Generation, Ext[1]], Index[Generation, Ext[2]],
  Index[Generation, Ext[3]], Index[Generation, Ext[4]]};

];

(* hack - manual correction of wrong sign generated by FeynRules for
vertices with 4 identical fermions and in uudd vertex *)

If [ SMEFT$Correct4FermionSign,

For[i=1, i < Length[FourLeptonVertices] + 1, i++, 
    If[ FourLeptonVertices[[i,1]] === {{lbar,1},{l,2},{lbar,3},{l,4}} ||
       (FourLeptonVertices[[i,1]] === {{vlbar,1},{vl,2},{vlbar,3},{vl,4}} &&
        ! SMEFT$MajoranaNeutrino),
        FourLeptonVertices[[i,2]] = FourLeptonVertices[[i,2]] /.
            TensDot[ia1_,ia2_][Index[Spin, Ext[1]], Index[Spin, Ext[4]]] ->
          - TensDot[ia1, ia2][Index[Spin, Ext[1]], Index[Spin, Ext[4]]]; 
    ];
];

For[i=1, i < Length[FourQuarkVertices] + 1, i++,
    If[ FourQuarkVertices[[i,1]] === {{uqbar,1},{uq,2},{uqbar,3},{uq,4}} ||
	FourQuarkVertices[[i,1]] === {{dqbar,1},{dq,2},{dqbar,3},{dq,4}} ||
	FourQuarkVertices[[i,1]] === {{dqbar,1},{dq,2},{uqbar,3},{uq,4}},
	FourQuarkVertices[[i,2]] = FourQuarkVertices[[i,2]] /.
	  TensDot[ia1_, ia2_][Index[Spin, Ext[1]], Index[Spin, Ext[4]]] ->
	- TensDot[ia1, ia2][Index[Spin, Ext[1]], Index[Spin, Ext[4]]];
	FourQuarkVertices[[i,2]] = FourQuarkVertices[[i,2]] /.
	  ProjM[Index[Spin, Ext[1]], Index[Spin, Ext[4]]] ->
	- ProjM[Index[Spin, Ext[1]], Index[Spin, Ext[4]]];
	FourQuarkVertices[[i,2]] = FourQuarkVertices[[i,2]] /.
	  ProjP[Index[Spin, Ext[1]], Index[Spin, Ext[4]]] ->
	- ProjP[Index[Spin, Ext[1]], Index[Spin, Ext[4]]];
    ];
  ];

,
(* if Correct4FermionSign = False and MajoranaNeutrino=True, revert
back manually evaluated proper sign to wrong one again ... *)

For[i=1, i < Length[FourLeptonVertices] + 1, i++, 
    If[ FourLeptonVertices[[i,1]] === {{vlbar,1},{vl,2},{vlbar,3},{vl,4}} &&
        SMEFT$MajoranaNeutrino,	
        FourLeptonVertices[[i,2]] = FourLeptonVertices[[i,2]] /.
	      TensDot[ia1_,ia2_][Index[Spin, Ext[1]], Index[Spin, Ext[4]]] ->
            - TensDot[ia1, ia2][Index[Spin, Ext[1]], Index[Spin, Ext[4]]]
    ];
];

];


find = Unique[fierz];

(* switch order of fields in Feynman rules for lvud *)

For[i=1, i < Length[TwoQuarkTwoLeptonVertices] + 1, i++,
    If[TwoQuarkTwoLeptonVertices[[i,1]] === {{dqbar,1},{l,2},{vlbar,3},{uq,4}},
       tmp = TwoQuarkTwoLeptonVertices[[i]][[2]] /. Ext[4] -> find;
       tmp = tmp /. Ext[2] -> Ext[4];
       TwoQuarkTwoLeptonVertices[[i]][[2]] = tmp /. find   -> Ext[2];
       TwoQuarkTwoLeptonVertices[[i]][[1]] = {{dqbar,1},{uq,2},{vlbar,3},{l,2}};
    ];

    If[TwoQuarkTwoLeptonVertices[[i,1]] === {{lbar,1},{dq,2},{uqbar,3},{vl,4}},
       tmp = TwoQuarkTwoLeptonVertices[[i]][[2]] /. Ext[4] -> find;
       tmp = tmp /. Ext[2] -> Ext[4];
       TwoQuarkTwoLeptonVertices[[i]][[2]] = tmp /. find    -> Ext[2];
       TwoQuarkTwoLeptonVertices[[i]][[1]] = {{lbar,1},{vl,2},{uqbar,3},{dq,2}};
    ];
];


];
(* end of FourFermionInteractions *)
];





BLViolating4FermionInteractions = Function[{},
(*  BL violatingfermion couplings *)
Module[{ oplist, i, j },

Print[Style["Calculating B and L violating 4-fermion vertices...",Bold]];

oplist = Intersection[ SMEFT$OperatorList, {"duq","qqu","qqq","duu"} ];
BLViolatingLagrangian = 0;
For [i=1, i<Length[oplist]+1, i++,
   BLViolatingLagrangian = BLViolatingLagrangian + Lam ToExpression[ "LQ"<>oplist[[i]] ];
];

BLViolatingLagrangian = ReplaceFlavorRotations[BLViolatingLagrangian] // Expand;
BLViolatingLagrangian = BLViolatingLagrangian // SMEFTExpOrder;
BLViolatingVertices = FeynmanRules[BLViolatingLagrangian];
BLViolatingVertices = BLViolatingVertices /.
       ToExpression[SMEFT$MB <> "qqu"][Index[Generation, Ext[1]], Index[Generation, Generation$1], i_, j_] -> 
       ToExpression[SMEFT$MB <> "qqu"][Index[Generation, Generation$1], Index[Generation, Ext[1]], i, j]; 
       
];
(* end of BLViolating4FermionInteractions *)
];




Neutrino4Vertex = Function[{},
(* if neutrinos are Majorana particles, 4-neutrino interaction must be
properly symmetrized, easier to do manually just for this one case *)

2 I Lam (ToExpression[SMEFT$MB<>"ll"][Index[Generation, Generation$1],
 Index[Generation, Generation$2], Index[Generation, Generation$3],
 Index[Generation, Generation$4]] ( Conjugate[ Ul[Index[Generation,
 Generation$1], Index[Generation, Ext[1]]] ] Ul[Index[Generation,
 Generation$2], Index[Generation, Ext[2]]] Conjugate[
 Ul[Index[Generation, Generation$3], Index[Generation, Ext[3]]] ]
 Ul[Index[Generation, Generation$4], Index[Generation, Ext[4]]]
 TensDot[Ga[Index[Lorentz, mu$1]], ProjM][Index[Spin, Ext[1]],
 Index[Spin, Ext[2]]] TensDot[Ga[Index[Lorentz, mu$1]],
 ProjM][Index[Spin, Ext[3]], Index[Spin, Ext[4]]] - Conjugate[
 Ul[Index[Generation, Generation$1], Index[Generation, Ext[2]]] ]
 Ul[Index[Generation, Generation$2], Index[Generation, Ext[1]]]
 Conjugate[ Ul[Index[Generation, Generation$3], Index[Generation,
 Ext[3]]] ] Ul[Index[Generation, Generation$4], Index[Generation,
 Ext[4]]] TensDot[Ga[Index[Lorentz, mu$1]], ProjP][Index[Spin,
 Ext[1]], Index[Spin, Ext[2]]] TensDot[Ga[Index[Lorentz, mu$1]],
 ProjM][Index[Spin, Ext[3]], Index[Spin, Ext[4]]] - Conjugate[
 Ul[Index[Generation, Generation$1], Index[Generation, Ext[1]]] ]
 Ul[Index[Generation, Generation$2], Index[Generation, Ext[2]]]
 Conjugate[ Ul[Index[Generation, Generation$3], Index[Generation,
 Ext[4]]] ] Ul[Index[Generation, Generation$4], Index[Generation,
 Ext[3]]] TensDot[Ga[Index[Lorentz, mu$1]], ProjM][Index[Spin,
 Ext[1]], Index[Spin, Ext[2]]] TensDot[Ga[Index[Lorentz, mu$1]],
 ProjP][Index[Spin, Ext[3]], Index[Spin, Ext[4]]] + Conjugate[
 Ul[Index[Generation, Generation$1], Index[Generation, Ext[2]]] ]
 Ul[Index[Generation, Generation$2], Index[Generation, Ext[1]]]
 Conjugate[ Ul[Index[Generation, Generation$3], Index[Generation,
 Ext[4]]] ] Ul[Index[Generation, Generation$4], Index[Generation,
 Ext[3]]] TensDot[Ga[Index[Lorentz, mu$1]], ProjP][Index[Spin,
 Ext[1]], Index[Spin, Ext[2]]] TensDot[Ga[Index[Lorentz, mu$1]],
 ProjP][Index[Spin, Ext[3]], Index[Spin, Ext[4]]] ) -
 ToExpression[SMEFT$MB<>"ll"][Index[Generation, Generation$1],
 Index[Generation, Generation$4], Index[Generation, Generation$3],
 Index[Generation, Generation$2]] ( Conjugate[ Ul[Index[Generation,
 Generation$1], Index[Generation, Ext[1]]] ] Ul[Index[Generation,
 Generation$4], Index[Generation, Ext[4]]] Conjugate[
 Ul[Index[Generation, Generation$3], Index[Generation, Ext[3]]] ]
 Ul[Index[Generation, Generation$2], Index[Generation, Ext[2]]]
 TensDot[Ga[Index[Lorentz, mu$1]], ProjM][Index[Spin, Ext[1]],
 Index[Spin, Ext[4]]] TensDot[Ga[Index[Lorentz, mu$1]],
 ProjM][Index[Spin, Ext[3]], Index[Spin, Ext[2]]] - Conjugate[
 Ul[Index[Generation, Generation$1], Index[Generation, Ext[4]]] ]
 Ul[Index[Generation, Generation$4], Index[Generation, Ext[1]]]
 Conjugate[ Ul[Index[Generation, Generation$3], Index[Generation,
 Ext[3]]] ] Ul[Index[Generation, Generation$2], Index[Generation,
 Ext[2]]] TensDot[Ga[Index[Lorentz, mu$1]], ProjP][Index[Spin,
 Ext[1]], Index[Spin, Ext[4]]] TensDot[Ga[Index[Lorentz, mu$1]],
 ProjM][Index[Spin, Ext[3]], Index[Spin, Ext[2]]] - Conjugate[
 Ul[Index[Generation, Generation$1], Index[Generation, Ext[1]]] ]
 Ul[Index[Generation, Generation$4], Index[Generation, Ext[4]]]
 Conjugate[ Ul[Index[Generation, Generation$3], Index[Generation,
 Ext[2]]] ] Ul[Index[Generation, Generation$2], Index[Generation,
 Ext[3]]] TensDot[Ga[Index[Lorentz, mu$1]], ProjM][Index[Spin,
 Ext[1]], Index[Spin, Ext[4]]] TensDot[Ga[Index[Lorentz, mu$1]],
 ProjP][Index[Spin, Ext[3]], Index[Spin, Ext[2]]] + Conjugate[
 Ul[Index[Generation, Generation$1], Index[Generation, Ext[4]]] ]
 Ul[Index[Generation, Generation$4], Index[Generation, Ext[1]]]
 Conjugate[ Ul[Index[Generation, Generation$3], Index[Generation,
 Ext[2]]] ] Ul[Index[Generation, Generation$2], Index[Generation,
 Ext[3]]] TensDot[Ga[Index[Lorentz, mu$1]], ProjP][Index[Spin,
 Ext[1]], Index[Spin, Ext[4]]] TensDot[Ga[Index[Lorentz, mu$1]],
 ProjP][Index[Spin, Ext[3]], Index[Spin, Ext[2]]] ) -
 ToExpression[SMEFT$MB<>"ll"][Index[Generation, Generation$1],
 Index[Generation, Generation$3], Index[Generation, Generation$2],
 Index[Generation, Generation$4]] ( Conjugate[ Ul[Index[Generation,
 Generation$1], Index[Generation, Ext[1]]] ] Ul[Index[Generation,
 Generation$3], Index[Generation, Ext[3]]] Conjugate[
 Ul[Index[Generation, Generation$2], Index[Generation, Ext[2]]] ]
 Ul[Index[Generation, Generation$4], Index[Generation, Ext[4]]]
 TensDot[Ga[Index[Lorentz, mu$1]], ProjM][Index[Spin, Ext[1]],
 Index[Spin, Ext[3]]] TensDot[Ga[Index[Lorentz, mu$1]],
 ProjM][Index[Spin, Ext[2]], Index[Spin, Ext[4]]] - Conjugate[
 Ul[Index[Generation, Generation$1], Index[Generation, Ext[3]]] ]
 Ul[Index[Generation, Generation$3], Index[Generation, Ext[1]]]
 Conjugate[ Ul[Index[Generation, Generation$2], Index[Generation,
 Ext[2]]] ] Ul[Index[Generation, Generation$4], Index[Generation,
 Ext[4]]] TensDot[Ga[Index[Lorentz, mu$1]], ProjP][Index[Spin,
 Ext[1]], Index[Spin, Ext[3]]] TensDot[Ga[Index[Lorentz, mu$1]],
 ProjM][Index[Spin, Ext[2]], Index[Spin, Ext[4]]] - Conjugate[
 Ul[Index[Generation, Generation$1], Index[Generation, Ext[1]]] ]
 Ul[Index[Generation, Generation$3], Index[Generation, Ext[3]]]
 Conjugate[ Ul[Index[Generation, Generation$2], Index[Generation,
 Ext[4]]] ] Ul[Index[Generation, Generation$4], Index[Generation,
 Ext[2]]] TensDot[Ga[Index[Lorentz, mu$1]], ProjM][Index[Spin,
 Ext[1]], Index[Spin, Ext[3]]] TensDot[Ga[Index[Lorentz, mu$1]],
 ProjP][Index[Spin, Ext[2]], Index[Spin, Ext[4]]] + Conjugate[
 Ul[Index[Generation, Generation$1], Index[Generation, Ext[3]]] ]
 Ul[Index[Generation, Generation$3], Index[Generation, Ext[1]]]
 Conjugate[ Ul[Index[Generation, Generation$2], Index[Generation,
 Ext[4]]] ] Ul[Index[Generation, Generation$4], Index[Generation,
 Ext[2]]] TensDot[Ga[Index[Lorentz, mu$1]], ProjP][Index[Spin,
 Ext[1]], Index[Spin, Ext[3]]] TensDot[Ga[Index[Lorentz, mu$1]],
 ProjP][Index[Spin, Ext[2]], Index[Spin, Ext[4]]] ) ) // Expand //
 Simplify

(* end of Neutrino4Vertex *)
]



