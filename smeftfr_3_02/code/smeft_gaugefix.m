(* SmeftFR v3.0 package *)
(* gauge fixing and Feynman rules for ghost terms (R_xi gauge only) *)

(*

   | W3 |          | Z |
   |    | = AZnorm |   |
   | B  |          | A |

*)

SMEFTLGhost = Function[{},
(* ghost Lagrangian in R_xi gauge *)
Module[{gvh,gvhbar,ghmat,ii,jj,k,aa,bb,cc,mu,xrot,xi,eta,aux,tmp,tmp1,tmp2},

xrot = Table[ AZnorm[ Index[SU2W,ii], Index[SU2W,jj] ],{ii,1,2},{jj,1,2}];

xi = Table[0,{ii,1,4},{jj,1,4}];
tmp = xrot.{{xiZ,0},{0,xiA}}.Transpose[xrot];
For[ii=1,ii<3,ii++,
  xi[[ii,ii]] = Wnorm^2 xiW;
  For[jj=1,jj<3,jj++, xi[[ii+2,jj+2]] = tmp[[ii,jj]];
  ];
];
xi = Expand[xi];

eta = Table[0,{ii,1,4},{jj,1,4}];
tmp = Inverse[xrot];
tmp = Transpose[tmp].tmp;
For[ii=1,ii<3,ii++,
  eta[[ii,ii]] = 1/Wnorm^2;
  For[jj=1,jj<3,jj++, eta[[ii+2,jj+2]] = tmp[[ii,jj]];
  ];		
];	
eta = Expand[eta];

gvh = {ghWi[1],ghWi[2],ghWi[3],ghB};
gvhbar = {ghWibar[1],ghWibar[2],ghWibar[3],ghBbar};
		
aux = gvhbar.eta.del[del[gvh,mu],mu];

aux = aux + gwnorm GW del[gvhbar,mu].eta.{{0,         -Wi[mu,3], Wi[mu,2],  0},
					  {Wi[mu,3],  0,         -Wi[mu,1], 0},
					  {-Wi[mu,2], Wi[mu,1],  0,         0},
					  {0,         0,         0,         0}}.gvh;

tmp = H/Hnorm + vev;
tmp1 = (GP + GPbar)/GPnorm/Sqrt[2];
tmp2 = I (GP - GPbar)/GPnorm/Sqrt[2];

ghmat = {
  {tmp,          G0/G0norm,   tmp1,         g1/gw tmp1    },
  {- G0/G0norm,  tmp,         tmp2,         g1/gw tmp2    },
  {- tmp1,       - tmp2,      tmp,          - g1/gw tmp   } G0norm^2,
  {g1/gw tmp1,   g1/gw tmp2,  - g1/gw tmp,  g1^2/gw^2 tmp } G0norm^2
};

ghmat = vev gw^2/4 Expand[xi.ghmat];
		
aux = aux + gvhbar.eta.ghmat.gvh // ExpandIndices;

aux = aux /. SMEFT$Identities;
		
(* QCD ghosts *)
aux = aux - del[ghGbar[aa],mu] del[ghG[aa],mu] + GS f[aa,bb,cc] del[ghGbar[aa],mu] G[mu,bb]  ghG[cc];

aux // ExpandIndices

]
(* end of SMEFTLGhost *)
];





SMEFTLGfix = Function[{},
(* gauge fixing terms in terms of PHYSICAL (mass eigenstate) A,Z,W,G0,G+ *)
Module[{aux,f1,f2,F,ixi,i,k,ii,jj,mu,nu,phfix,pzfix,pwfix,pzf1,pzf2,pzf3,xrot,xi,VEV,tmp},
	
VEV = {0,vev/Sqrt[2]};

xrot = {{gwnorm,0},{0,g1norm}}.Table[ AZnorm[ Index[SU2W,ii], Index[SU2W,jj] ],{ii,1,2},{jj,1,2}] // Expand;

xi = Table[0,{ii,1,4},{jj,1,4}];
tmp = xrot.{{xiZ,0},{0,xiA}}.Transpose[xrot];
For[ii=1,ii<3,ii++,
  xi[[ii,ii]] = xiW;
  For[jj=1,jj<3,jj++, xi[[ii+2,jj+2]] = tmp[[ii,jj]];
  ];		
];	
xi = xi // Simplify;

f2 = Table[0,{ii,1,4}];
ixi = Inverse[xi] // Simplify;

f1 = {ExpandIndices[gwnorm del[Wi[mu,1],mu]],
      ExpandIndices[gwnorm del[Wi[mu,2],mu]],
      ExpandIndices[gwnorm del[Wi[mu,3],mu]],
      ExpandIndices[g1norm del[B[mu],mu]]};

For[i=1,i<4,i++,
    For[ii=1,ii<3,ii++,
	For[jj=1,jj<3,jj++,
	    f2[[i]] = f2[[i]] + GW Ta[i,ii,jj] (HC[Phi[ii]] VEV[[jj]] - VEV[[ii]] Phi[jj]);
	];
    ];
];
f2[[3]] = G0norm^2 f2[[3]] ;

For[ii=1,ii<3,ii++, f2[[4]] = f2[[4]] + (HC[Phi[ii]] VEV[[ii]] - Phi[ii] VEV[[ii]]) ];
f2[[4]] = G1/2 G0norm^2 f2[[4]];

f2 = ExpandIndices[f2, FlavorExpand->SU2W] // Simplify;

F = f1 - xi.f2;

aux = - F.ixi.(F /. mu -> nu)/2  // Expand;

phfix = aux - (aux /. A[___]->0);
aux = aux - phfix // Simplify;
pzfix = aux - (aux /. Z[___]->0 /. G0->0);
pwfix = aux - pzfix;

phfix = phfix /. SMEFT$Identities // Simplify;
pzfix = pzfix /. SMEFT$Identities // Simplify;
pwfix = pwfix /. SMEFT$Identities // Simplify;

phfix + pzfix + pwfix	

]
(* end of SMEFTLGfix *)
];


GhostInteractions = Function[{},
(* EW and QCD ghost couplings *)

Print[Style["Calculating ghost vertices...",Bold]];
(* run it to check the proper form of gauge-fixing terms
   test = SMEFTLGfix[]; *)

GhostLagrangian = (SMEFTLGhost[] /. SMEFT$Identities) // OptimizeIndex// Simplify;
GhostVertices = FeynmanRules[GhostLagrangian] // Simplify;

];

