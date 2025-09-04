(* SmeftFR v3.0 package *)
(* routines for manipulating/updating SmeftFR parameters and internal
variables *)



CheckClashingNames = Function[{},
(*  test if user defined parameter names does not overlap with SmeftFR default variables *)
Module[{user, smeft, tmp},

user = ToString /@ SM$InputParameters[[All,1]];
If[ ! DuplicateFreeQ[user],
  Print[Style["Double defined user parameter names, check header of smeft_input_scheme.m!\n"],Bold];
  Abort[];
];

(* SmeftFR parameters *)
tmp = ToExpression[ReadString[ FileNameJoin[{SMEFT$Path, "definitions", "smeft_par_head_WB.fr"}] ] <> "{xxx}}"];
tmp = Drop[tmp[[All, 1]],-1];
smeft = ToString /@ tmp;
tmp = ToExpression[ReadString[ FileNameJoin[{SMEFT$Path, "definitions", "smeft_par_head_MB.fr"}] ] <> "{xxx}}"];
tmp = Drop[tmp[[All, 1]],-1];
smeft = Join[smeft, ToString /@ tmp];
Get[ FileNameJoin[{SMEFT$Path, "definitions", "smeft_par_SM.fr"}] ];
smeft = Join[smeft, ToString /@ SMEFT$SMParameters[[All,1]] ];

(* SmeftFR masses *)
Get[ FileNameJoin[{SMEFT$Path, "definitions", "smeft_fields_WB.fr"}] ];
tmp = Flatten[Mass /. Values[Rule @@@ M$ClassesDescription] ];
tmp = Complement[ ToString /@ tmp,  {"0","Mass","Internal", "External"} ];
smeft = Join[smeft, tmp];
Get[ FileNameJoin[{SMEFT$Path, "definitions", "smeft_fields_MB.fr"}] ];
tmp = Flatten[Mass /. Values[Rule @@@ M$ClassesDescription] ];
tmp = Complement[ ToString /@ tmp,  {"0","Mass","Internal","External"} ];
smeft = Join[smeft, tmp];

(* SmeftFR operators *)
tmp = (SMEFT$WB <> #) & /@ SMEFT$AllOperators;
smeft = Join[smeft, tmp];
tmp = (SMEFT$MB <> #) & /@ SMEFT$AllOperators;
smeft = DeleteDuplicates[Join[smeft, tmp]];

tmp = Intersection[user, smeft];
If [ Length[tmp] > 0,
   Print[Style["User defined parameter names in the header of smeft_input_scheme.m clash with standard SmeftFR parameters. Change names of the following parameters: ",Bold], tmp];
   Abort[];
];

]
(* end of CheckClashingNames *)			
]




UpdateSMParameters = Function[{},
(* evaluate numerical values of parameters in SM Lagrangian *)

Module[{num, stdl, mlist, tmp, sminput, smnum, smpar, smnames,
smvalues, x, y},

num = Rule @@@ Partition[ Riffle[ SM$InputParameters[[All, 1]],
                         Value /. SM$InputParameters[[All, 2]]], 2];

(* below GS is commented out as MG5 requires it to be defined as
   "internal" and expressed in terms of aS *)

stdl = {"G1", "GW", (* "GS",*) "vev", "hlambda", "MZ", "MW" , "MH",
"MQT", "MQC", "MQU", "MQB", "MQS", "MQD", "MLT", "MLM", "MLE", "MVT",
"MVM", "MVE"};

mlist = {"mz", "mw" , "mh", "mqt", "mqc", "mqu", "mqb", "mqs", "mqd",
  "mlt", "mlm", "mle", "mvt", "mvm", "mve"};

tmp = "UserInput$" <> # & /@ stdl;
tmp = StringDrop[#,StringLength["UserInput$"]] & /@ Complement[ tmp, Intersection[ToString /@ (ToExpression /@ tmp), tmp] ];
sminput = Rule @@@ Partition[ Riffle[ ToExpression /@ tmp, ToExpression["UserInput$" <> #] & /@ tmp], 2];

smnum = ToExpression["SMnum$" <> ToString[#]] -> ToExpression["UserInput$" <> ToString[#]] & /@ stdl;
smnum = smnum //. Join[num, SMEFT$WCValues, sminput];

If[ ! NumberQ[Values[#]],
  Print[Style["Numerical value of " <> StringDrop[ToString[Keys[#]],6] <>
              " not fully initialized, check completeness of the input scheme!", Bold, Red]];
  Abort[];
  ] & /@ smnum;

smnum = Association[smnum /. Complex[x_, 0.] -> x];

If[ ! NumberQ[#],
  Print[Style["Numerical values of CKM matrix elements  not fully initialized, check completeness of the input scheme!",Bold,Red]]; Abort[]; ] & /@ Flatten[UserInput$CKM];

If[ ! NumberQ[#],
  Print[Style["Numerical values of PMNS matrix elements  not fully initialized, check completeness of the input scheme!",Bold,Red]]; Abort[]; ] & /@ Flatten[UserInput$PMNS];

smnum =  Append[smnum,  { SMnum$Kq -> Flatten[ Table[ Kq[i,j] -> UserInput$CKM[[i,j]], {i,1,3},  {j,1,3} ] ],
                          SMnum$Ul -> Flatten[ Table[ Ul[i,j] -> UserInput$PMNS[[i,j]], {i,1,3}, {j,1,3} ] ] } ];

(* update SM parameters *)
Get[ FileNameJoin[{SMEFT$Path, "definitions", "smeft_par_SM.fr"}] ];
smpar = Association[SMEFT$SMParameters /. Equal -> Rule];
smnames = Intersection[ToString /@ Keys[smpar], Join[stdl,{"Kq","Ul","AZnorm"}]];
smvalues = smpar[ToExpression[#]] /. (Value -> x_) -> (Value -> smnum[ToExpression["SMnum$" <> #]]) & /@ smnames;
smvalues = Association[Rule @@@ Partition[ Riffle[ ToExpression /@ smnames, smvalues], 2]];
tmp = (ToExpression[#] == x_) -> (ToExpression[#] == smvalues[ToExpression[#]]) & /@ smnames;
SMEFT$SMParameters = SMEFT$SMParameters /. tmp;

(* update SM particle masses *)
smnames = Intersection[ToString /@ Keys[smpar], mlist];
smvalues = smpar[ToExpression[#]] /. (Value -> x_) ->
       (Value -> smnum[ ToExpression[ "SMnum$" <> ToUpperCase[#] ] ] ) & /@ smnames;
smvalues = Association[Rule @@@ Partition[ Riffle[ ToExpression /@ smnames, smvalues], 2]];
tmp = (ToExpression[#] == x_) -> (ToExpression[#] == smvalues[ToExpression[#]]) & /@ smnames;
SMEFT$SMParameters = SMEFT$SMParameters /. tmp;

SMEFT$SMParametersReal = SMEFT$SMParameters /. Complex[x_, y_] -> x;

]
(* end of UpdateSMParameters *)			
];






UpdateNormalizationConstants = Function[{},
(* test if user defined parameter names does not overlap with SmeftFR
variables *)

Module[{num, stdl, nlist, tmp, sminput, smnum, smpar, smnames,
smvalues, x},

num = Rule @@@ Partition[ Riffle[ SM$InputParameters[[All, 1]],
                         Value /. SM$InputParameters[[All, 2]]], 2];

stdl = {"G1", "GW", "GS", "vev", "hlambda", "MZ", "MW" , "MH",
  "MQT", "MQC", "MQU", "MQB", "MQS", "MQD", "MLT", "MLM", "MLE",
  "MVT", "MVM", "MVE"};

nlist = {"gsnorm", "g1norm", "gwnorm", "Hnorm", "G0norm", "GPnorm", "Wnorm", "Gnorm", "AZnorm11",  "AZnorm21",  "AZnorm12",  "AZnorm22"};

SMEFT$AZnorm11 = SMEFT$AZnorm[[1, 1]];
SMEFT$AZnorm12 = SMEFT$AZnorm[[1, 2]];
SMEFT$AZnorm21 = SMEFT$AZnorm[[2, 1]];
SMEFT$AZnorm22 = SMEFT$AZnorm[[2, 2]];

tmp = "UserInput$" <> # & /@ stdl;
tmp = StringDrop[#,StringLength["UserInput$"]] & /@ Complement[ tmp, Intersection[ToString /@ (ToExpression /@ tmp), tmp] ];
sminput = Rule @@@ Partition[ Riffle[ ToExpression /@ tmp, ToExpression["UserInput$" <> #] & /@ tmp], 2];

(* numerical values of normalization constants *)
smnum = ToExpression["SMnum$" <> ToString[#]] -> SMEFTExpOrder[ ToExpression["SMEFT$" <> #] /. sminput ] & /@ nlist //. Join[num, SMEFT$WCValues, sminput];

smnum = Association[smnum];

tmp = Table[ SMEFTExpOrder[ SMEFT$AZnorm /. sminput ][[i,j]], {i,1,2},  {j,1,2} ] //. Join[num, SMEFT$WCValues, sminput];
smnum  =  Append[smnum,  { SMnum$AZnorm -> Flatten[ Table[ AZnorm[i,j] -> tmp[[i,j]], {i,1,2},  {j,1,2} ] ]} ];

(* update normalization constants in the list *)
smpar = Association[SMEFT$SMParameters /. Equal -> Rule];
smnames = Intersection[ToString /@ Keys[smpar], Join[nlist,{"AZnorm"}]];
smvalues = smpar[ToExpression[#]] /. (Value -> x_) -> (Value -> smnum[ToExpression["SMnum$" <> #]]) & /@ smnames;
smvalues = Association[Rule @@@ Partition[ Riffle[ ToExpression /@ smnames, smvalues], 2]];
tmp = (ToExpression[#] == x_) -> (ToExpression[#] == smvalues[ToExpression[#]]) & /@ smnames;
SMEFT$SMParameters = SMEFT$SMParameters /. tmp;
SMEFT$SMParametersReal = SMEFT$SMParametersReal /. tmp;

]
(* end of UpdateNormalizationConstants *)			
];




NormalizationConstantsExpansion = Function[{input},
(* reexpress normalization constants in terms of user-defined input quantities *)
Module[{norm, const, azind, aznorm, clist, nlist, alist, name, a, b, c},

If [ ! SMEFT$UserInputInitialization && ! SMEFT$InputListStatus,
     Print["User input scheme not initialized - call SMEFTInputScheme[] first (modifying it if necessary)"];
     Abort[];
   ];
 
const = {"GW", "G1", "GS", "hlambda", "vev", "MZ", "MW", "MH"};
norm  = {"g1norm", "gwnorm", "gsnorm", "Wnorm", "Gnorm", "Hnorm", "G0norm", "GPnorm"};
azind = {{1,1},{1,2},{2,1},{2,2}};

aznorm = "AZnorm[Index[SU2W, " <> ToString[ #[[1]] ] <>
             "], Index[SU2W, " <> ToString[ #[[2]] ]  <> "] ]" & /@ azind;


SMEFT$InputList = (ToExpression[#] -> ToExpression["Dim4" <> #] +
                                      ToExpression["Dim6" <> #] Lam +
                                      ToExpression["Dim8" <> #] Lam^2) & /@ Join[const,aznorm,norm];

SMEFT$InputList = SMEFT$InputList /. (ToExpression["Dim4" <> #] -> 1 & /@ norm);

Switch[ ToUpperCase[input],
  "USER",  name = "UserInput$";
           SMEFT$InputSub = User$InputSub,
  "SMEFT", name = "SMEFT$";
           SMEFT$InputSub = Default$InputSub,
  "NONE",  name = "NONE";
	   SMEFT$InputSub = {};
	   SMEFT$InputList = {},
   _,      Print["Argument of NormalizationConstantsExpansion should be \"user\", \"smeft\" or \"none\""];
           Abort[];
];

       
If [ ! SMEFT$InputListStatus && name != "NONE",
  clist = ToExpression[#] -> ToExpression[name <> # ] & /@ const;
  nlist = ToExpression[#] -> SMEFTExpOrder[ ToExpression["SMEFT$" <> #] /. clist ] & /@ norm;
  alist = SMEFTExpOrder[ SMEFT$AZnorm[[ ToExpression[#][[1]], ToExpression[#][[2]] ]] /. clist ] & /@ azind;
  alist = Rule @@@ Partition[ Riffle[ToExpression /@ aznorm, alist ], 2 ];   

  nlist = Join[clist, nlist, alist];
       
  clist = Keys[nlist];
  alist = Values[nlist] //. nlist;
  nlist = Rule @@@ Partition[ Riffle[clist, alist], 2];

  SMEFT$InputSub = (ToExpression["Dim4" <> ToString[Keys[#]]] ->
                   (Expand[Values[#]] /. Lam -> 0 // Simplify)) & /@ nlist;
  SMEFT$InputSub = Complement[ SMEFT$InputSub, (ToExpression["Dim4" <> #] -> 1 & /@ norm) ];

  clist          = (ToExpression["Dim6" <> ToString[Keys[#]]] ->
                   (Coefficient[Expand[Values[#]],Lam] // Simplify)) & /@ nlist;

  alist          = (ToExpression["Dim8" <> ToString[Keys[#]]] ->
                   (Coefficient[Expand[Values[#]],Lam^2] // Simplify)) & /@ nlist;

  SMEFT$InputSub = Join[ SMEFT$InputSub, clist, alist ] // FunctionExpand; 
  SMEFT$InputSub = SMEFT$InputSub //. (a_ b_)^c_ -> a^c b^c;
  SMEFT$InputSub = SMEFT$InputSub //. (a_/b_)^c_ -> a^c/b^c;
  SMEFT$InputSub = SMEFT$InputSub //. (1/a_)^b_ -> 1/a^b;
  SMEFT$InputSub = SMEFT$InputSub /. (UserInput$MW^a_)^b_ -> UserInput$MW^(a b);
  SMEFT$InputSub = SMEFT$InputSub /. (UserInput$MZ^a_)^b_ -> UserInput$MZ^(a b);
  SMEFT$InputSub = SMEFT$InputSub // Simplify;
  SMEFT$InputSub = SMEFT$InputSub //. (a_ b_)^c_ -> a^c b^c;

(* MG5 ugly hack: GS should be always expressed in terms of aS
parameter; aS is evaluated in UpdateSpecialParameters using values of
input parameters... *) 
  If[ name == "UserInput$",
    SMEFT$InputSub = SMEFT$InputSub /. (Dim4GS -> a_) -> (Dim4GS -> 2 Sqrt[Pi] Sqrt[aS]) /.
                                       (Dim6GS -> a_) -> (Dim6GS -> 0) /.
                                       (Dim8GS -> a_) -> (Dim8GS -> 0);
  ];       
];
 
]
(* end of NormalizationConstantsExpansion *)
]



OutputParametersExpansion = Function[{expr},
(* expand Feynman rules substituting SM shift constants with explicit 1/Lambda^n terms to a given order *)
Module[{tmp,tmp1},

Switch[ SMEFT$ExpansionOrder,
  0, tmp = (expr /. Lam->0) //. (SMEFT$InputList /. Lam->0) // Expand,
  1, tmp = OutputParametersExpansionDim6[expr],
  2, tmp = OutputParametersExpansionDim8[expr],
  _, Print["SMEFT$ExpansionOrder in OutputParametersExpansion should be 0, 1 or 2, it is ", SMEFT$ExpansionOrder];
     Abort[];
];

(* keep only terms linear in neutrino mass *)
tmp = tmp /. Lam fmv[__] -> 0 /. Lam^2 fmv[__] -> 0;

(* replace generic expanded Z_X by expressions in terms of WCs *)       
tmp = tmp /. SMEFT$InputSub;  

(*simplify the result *)       
If [(Head[expr] != List) || (Length[expr] == 0),
  tmp = CombineCommonPowers[tmp, Lam, SMEFT$ExpansionOrder];
,
  If[ tmp[[1,1,2]] === 1,
     tmp[[2]] = CombineCommonPowers[tmp[[2]], Lam, SMEFT$ExpansionOrder];
  ,
    tmp1 = CombineCommonPowers[#, Lam, SMEFT$ExpansionOrder] & /@ tmp[[All, 2]];
    tmp = Partition[Riffle[tmp[[All,1]], tmp1],2];
  ];
];

tmp       
       
]
(* end of OutputParametersExpansion *)
]
				     


OutputParametersExpansionDim8 = Function[{expr},
(* expand Feynman rules substituting SM shift constants with explicit max 1/Lambda^4 terms *)
Module[{tmp,tmp1,tmp2},

tmp  = (expr /. Lam->0) /. SMEFT$InputList // SMEFTExpOrder;

tmp1 = Coefficient[Expand[expr],Lam] /. (SMEFT$InputList /. Lam^2->0) // Expand; 
tmp1 = tmp1 /. Lam^k_ -> If[k>1,0,Lam^k];
tmp1 = Normal[Series[tmp1, {Lam,0,1}]];

tmp2 = Coefficient[Expand[expr],Lam^2] /. (SMEFT$InputList /. Lam->0) // Expand;
tmp2 = tmp2 /. Lam->0;

tmp + Lam tmp1 + Lam^2 tmp2 // Expand

]
(* end of OutputParametersExpansionDim8 *)			
]




OutputParametersExpansionDim6 = Function[{expr},
(* expand Feynman rules substituting SM shift constants with explicit max 1/Lambda^2 terms *)
Module[{tmp,tmp1},

tmp  = (expr /. Lam->0) //. (SMEFT$InputList /. Lam^2->0) // SMEFTExpOrder;
tmp  = tmp /. Lam^k_ -> If[k>1,0,Lam^k]; 

tmp1 = Coefficient[Expand[expr],Lam] //. (SMEFT$InputList /. Lam->0) // Expand;  

tmp + Lam tmp1 // Expand

]
(* end of OutputParametersExpansionDim6 *)			
]




UpdateSpecialParameters = Function[{},
(* evaluate analytical expressions and numerical values of "special"
parameters which names are fixed (hard-coded) in external programs *)

Block[{vev, Gf}, Module[{num, special, g1, gw, eps, tmp, sinput, snum,
spar, snames, svalues, mdiff, simplify},

num = Rule @@@ Partition[ Riffle[ SM$InputParameters[[All, 1]],
                         Value /. SM$InputParameters[[All, 2]]], 2];

(* currently the list below contains parameters used by Madgraph 5;
   add more if necessary

Entries for parameters listed below has to be also added in the file
definitions/smeft_par_SM.fr, their values are only updated here! *)

special = {"aS", "aEWM1"};

UserInput$aS = UserInput$GS^2/4/Pi // SMEFTExpOrder;       

gw  = UserInput$GW;
g1  = UserInput$G1;
eps = Lam vev^2/SMEFT$gwnorm/SMEFT$g1norm (MBScalarList["phiWB"] +
                               Lam vev^2/2 MBScalarList["WBphi4n1"]) /.
                                   vev -> UserInput$vev// SMEFTExpOrder;

UserInput$aEWM1 = g1^2 gw^2/(g1^2 + 2 eps g1 gw + gw^2)/4/Pi // SMEFTExpOrder;
  
(* some simplifications, could be scheme dependent... *)
simplify = Module[{tmp,a,b,k},
  tmp = # /. UserInput$MZ -> Sqrt[mdiff^2 + UserInput$MW^2] // Expand;
  tmp = Simplify[tmp, Assumptions -> Gf > 0 && mdiff > 0 && UserInput$MW > 0];
  tmp = tmp //. (a_ b_)^k_ -> a^k b^k;
  tmp = tmp /. mdiff -> Sqrt[UserInput$MZ^2 - UserInput$MW^2];
  tmp = Simplify[tmp, Assumptions -> Gf > 0 && UserInput$MZ > 0 && UserInput$MW > 0 && UserInput$MZ > UserInput$MW];
  tmp = CombineCommonPowers[tmp, Lam, SMEFT$ExpansionOrder] //. (a_ b_)^k_ -> a^k b^k;
  tmp /. (1/a_)^k_ -> 1/a^k	  
] & ;

UserInput$aS    = simplify[UserInput$aS];
UserInput$aEWM1 = simplify[UserInput$aEWM1];
			
(* update numerical values of variables in the list *)
tmp = "UserInput$" <> # & /@ special;
tmp = StringDrop[#,StringLength["UserInput$"]] & /@ Complement[ tmp, Intersection[ToString /@ (ToExpression /@ tmp), tmp] ];
sinput = Rule @@@ Partition[ Riffle[ ToExpression /@ tmp, ToExpression["UserInput$" <> #] & /@ tmp], 2]; 
      
snum = ToExpression["SMnum$" <> ToString[#]] -> ToExpression["UserInput$" <> ToString[#]] & /@ special;
snum = snum //. Join[num, SMEFT$WCValues, sinput];

If[ ! NumberQ[Values[#]],
  Print[Style["Numerical value of " <> StringDrop[ToString[Keys[#]],6] <>
              " not fully initialized, check completeness of the input scheme!", Bold, Red]];
  Abort[];
  ] & /@ snum;

snum = Association[snum /. Complex[x_, 0.] -> x];

spar = Association[SMEFT$SMParameters /. Equal -> Rule];
snames = Intersection[ToString /@ Keys[spar], special];
svalues = spar[ToExpression[#]] /. (Value -> x_) -> (Value -> snum[ToExpression["SMnum$" <> #]]) & /@ snames;
svalues = Association[ Rule @@@ Partition[ Riffle[ ToExpression /@ snames, svalues], 2] ];
tmp = (ToExpression[#] == x_) -> (ToExpression[#] == svalues[ToExpression[#]]) & /@ snames;    
SMEFT$SMParameters = SMEFT$SMParameters /. tmp;
SMEFT$SMParametersReal = SMEFT$SMParametersReal /. tmp;


]]
(* end of UpdateSpecialParameters *)			
];



