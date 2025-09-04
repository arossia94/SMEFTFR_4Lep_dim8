(* SmeftFR v3.0 package *)
(* latex output generator *)

(* latex output auxiliary variables *)

FieldType = {H,G0,GP,GPbar,G,Z,A,W,Wbar,
	     l,vl,dq,uq,lbar,vlbar,dqbar,uqbar,
             ghWp,ghWm,ghZ,ghA,ghG,
	     ghWpbar,ghWmbar,ghZbar,ghAbar,ghGbar,
             lc,vlc,dqc,uqc,lcbar,vlcbar,dqcbar,uqcbar};

InFieldName = {"h","G0","GP","GP","g","Z","A","W","W",
	     "e","v","dq","uq","ec","vc","dqc","uqc",
             "etaP","etaM","etaZ","etaA","etaG",
             "etabarP","etabarM","etabarZ","etabarA","etabarG",
             "ec","v","dq","uq","ec","v","dq","uq"};

OutFieldName = {"h","G0","GP","GP","g","Z","A","W","W",
	     "ec","vc","dqc","uqc","e","v","dq","uq",
             "etaP","etaM","etaZ","etaA","etaG",
             "etabarP","etabarM","etabarZ","etabarA","etabarG",
             "ec","v","dq","uq","ec","v","dq","uq"};

InFieldName  = Rule @@@ Partition[Riffle[FieldType, InFieldName],2];
OutFieldName = Rule @@@ Partition[Riffle[FieldType, OutFieldName],2];

FieldTex = {"h","G^0","G^+","G^-","g","Z^0","A^0","W^+","W^-",
            "e","\\nu","d","u"," e","\\nu"," d"," u",
           "\\eta^+","\\eta^-","\\eta_Z","\\eta_A","\\eta_G",
	   "\\bar\\eta^+","\\bar\\eta^-","\\bar\\eta_Z","\\bar\\eta_A","\\bar\\eta_G",
           "(e^C)","(\\nu^C)","(d^C)","(u^C)","(e^C)","(\\nu^C)","(d^C)","(u^C)"};

LineType =  {1,1,2,2,3,4,4,7,7,
	     6,6,6,6,5,5,5,5,
             8,8,8,8,8,
	     9,9,9,9,9,
             6,6,6,6,5,5,5,5};


(* auxiliary functions for generation of diagrams with latex version
of Feynman rules *)


Dim6VertexName = Function[{plist},
(* find vertex name tag *)
Module[{tmp,f,fc},
  f = {"v","e","uq","dq"};
  fc = {"vc","ec","uqc","dqc"};
  tmp = plist[[All,1]];				 
  tmp[[1]] = tmp[[1]] /. OutFieldName;
  tmp[[2]] = tmp[[2]] /. InFieldName; 
  tmp[[3]] = tmp[[3]] /. OutFieldName;
  If [Length[tmp] > 3, tmp[[4]] = tmp[[4]] /. InFieldName];
  tmp = tmp /. OutFieldName; 
  If [Length[tmp] == 4 && MemberQ[fc,tmp[[1]]] && MemberQ[f,tmp[[2]]] &&
      MemberQ[fc,tmp[[3]]] && MemberQ[f,tmp[[4]]],
    tmp[[1]] = StringDrop[tmp[[1]],-1];
    tmp[[3]] = StringDrop[tmp[[3]],-1];
    tmp[[2]] = tmp[[2]] <> "c";
    tmp[[4]] = tmp[[4]] <> "c";
  ];
  StringJoin[ Sort[ tmp ] ]
]

];


FieldNameAndIndex = Function[{pos,i,DiagType},
Module[{flav,ind,is},

flav  = {0,0,0,0,3,4,4,4,4,1,1,2,2,1,1,2,2,0,0,0,0,5,0,0,0,0,5,1,1,2,2,1,1,2,2};
is = ToString[i];

ind = Switch[flav[[pos]],
	     0, "",
	     1, "^{f_" ~~ is ~~ "}" ~~ If[DiagType == 2, "_{s_" ~~ is ~~ "}",""],
	     2, "^{f_" ~~ is ~~ "}" ~~ Switch[ DiagType,
					       1, "_{m_" ~~ is ~~ "}",
					       2, "_{m_" ~~ is ~~ " s_" ~~ is ~~ "}",
					       _, "" ],
	     3, "^{a_" ~~ is ~~ "}_{\\mu_" ~~ is ~~ "}",
	     4, "_{\\mu_" ~~ is ~~ "}",
	     5, "^{a_" ~~ is ~~ "}"
];

"$" ~~ FieldTex[[pos]] ~~ ind ~~ "$"

]
(*end of FieldNameAndIndex *)
]



IsolateFactors = Function[{expr},
Module[{flist,f1,f2,i},

f1 = f2 = 1;
flist = FactorList[Simplify[expr, TimeConstraint->10]] ;
For[i=1, i < Length[flist] + 1, i++,
    If[ ! MemberQ[ SMEFT$Dim6BosonOperators, flist[[i,1]] ] &&
	  MemberQ[ {Power, Complex, Symbol, Integer}, Head[ flist[[i,1]] ] ],
	f1 = f1 flist[[i,1]]^flist[[i,2]];
	,
	f2 = f2 flist[[i,1]]^flist[[i,2]];
    ];
];

{f1,f2}

]
(* end of IsolateFactors *)
]



ToLaTeXForm = Function[{expr},
Module[{sm,nsm,op,coeff,lorentz,res,num1,num2,den,i,l,a,AAAA,AAAB,AAA},

l = Length[SMEFT$Dim6NullList];
a = $ModuleNumber - 1;

sm = expr /. SMEFT$Dim6NullList;
nsm = expr - sm // Expand;

coeff = {};
lorentz = {};

sm = Together[sm // Simplify];
If[ sm =!= 0,
    den = Denominator[sm];
    {num1,num2} = IsolateFactors[Numerator[sm]];
(*   num1 = num1 /. G1 - GW -> (G1^2 - GW^2)/(G1 + GW) /. GW - G1 -> - (G1^2 - GW^2)/(G1 + GW);*)
    coeff = AppendTo[coeff, Simplify[num1, TimeConstraint->10]/den];
    lorentz = AppendTo[lorentz, Simplify[num2, TimeConstraint->10]];
];

For[i=1, i < l+1, i++,
  op = nsm /. Drop[SMEFT$Dim6NullList,{i}];
  op = Together[ Simplify[op, TimeConstraint->5] ];
  If[ op =!= 0,
    den = Denominator[op];
    {num1,num2} = IsolateFactors[Numerator[op]];
(*   num1 = num1 /. G1 - GW -> (G1^2 - GW^2)/(G1 + GW) /. GW - G1 -> - (G1^2 - GW^2)/(G1 + GW);*)
    coeff = AppendTo[coeff, Simplify[num1, TimeConstraint->10]/den];
    lorentz = AppendTo[lorentz, Simplify[num2, TimeConstraint->10]];
  ];
];

coeff = coeff /. G1-> AAAA /. GW->AAAB /. mvdiag[i_] -> vmass[i];
coeff = ToString[ TeXForm[ # ] ] & /@ coeff;
(* remove blank spaces *)
coeff = StringReplace[coeff, (StartOfString ~~ Whitespace) | (Whitespace ~~ EndOfString) :> ""];
For[ i=1, i < Length[coeff]+1, i++,
  If [ coeff[[i]] == "1",  coeff[[i]] = " " ];
  If [ coeff[[i]] == "-1" || coeff[[i]] == "- 1",  coeff[[i]] = "-" ];
];
(* set nicely formatted first sign in the line *)
coeff = If[ StringTake[#,1] == "-", " - " <> StringDrop[#, 1], " + " <> # ] & /@ coeff;

lorentz = lorentz /. G1-> AAAA /. GW->AAAB /. mvdiag[i_] -> vmass[i];
For[ i=1, i<Length[lorentz]+1, i++, If [ Head[ lorentz[[i]] ] === Plus, lorentz[[i]] = AAA lorentz[[i]] ] ];
lorentz = ToString[ TeXForm[ # ] ] & /@ lorentz;
For[ i=1, i < Length[lorentz]+1, i++,
  If [ lorentz[[i]] == "1",  lorentz[[i]] = " " ];
  If [ lorentz[[i]] == "-1" || lorentz[[i]] == "- 1",  lorentz[[i]] = "-" ];
];

res = Riffle[coeff,lorentz];
res = StringReplace[res, "AAAA$\\$$"<>ToString[a] -> "AAAA"];
res = StringReplace[res, "AAAB$\\$$"<>ToString[a] -> "AAAB"];
res = StringReplace[res, "AAA$\\$$" <>ToString[a] -> ""];

res

]
(* end of ToLaTexForm *)
]




ToLaTeX = Function[{expr,DiagType,verbose},
Block[{SIG},Module[{tmp,res,x,y,v,z,a,b,c,i,s1,s2,s3,s4,sprod},

tmp = expr /. Lam -> 1;
tmp = tmp  /. TensDot[Ga[x_], SlashedP[y_], z_][s1_, s2_] ->
              FV[y, x] z[s1,s2] - I SIG[x, \[Nu], z, s1, s2] FV[y, \[Nu]]  /.
              TensDot[SlashedP[y_], Ga[x_], z_][s1_, s2_] ->
              FV[y, x] z[s1,s2] + I SIG[x, \[Nu], z, s1, s2] FV[y, \[Nu]]  // Expand;

tmp = tmp /. TensDot[Ga[x_], Ga[y_], z_][s1_, s2_] -> ME[x, y] z[s1,s2] - I SIG[x, y, z, s1, s2] // Expand;

tmp = tmp /. SIG[Index[Lorentz, mu$2], Index[Lorentz, mu$1], x_, s1_, s2_] ->
           - SIG[Index[Lorentz, mu$1], Index[Lorentz, mu$2], x, s1, s2];

tmp = tmp /. SIG[Index[Lorentz, Ext[a_]], Index[Lorentz, Ext[b_]], x_, y_, z_] ->
    If[a <= b,  SIG[Index[Lorentz, Ext[a]], Index[Lorentz, Ext[b]], x, y, z],
              - SIG[Index[Lorentz, Ext[b]], Index[Lorentz, Ext[a]], x, y, z] ];
    
tmp = tmp /. Eps[x_, y_, v_, z_] FV[a_, y_] FV[b_, v_] FV[c_, z_] ->
    Eps[x, y, v, z] Switch[{a > b, a > c, b > c},
	   {False, False, False} ,  FV[a, y] FV[b, v] FV[c, z],
	   {False, False, True}, -  FV[a, y] FV[b, z] FV[c, v],
	   {True, False, False}, -  FV[a, v] FV[b, y] FV[c, z],
	   {True, True, True},   -  FV[a, z] FV[b, v] FV[c, y],
	   {True, True, False},     FV[a, z] FV[b, y] FV[c, v],
	   {False, True, True},     FV[a, v] FV[b, z] FV[c, y]];

tmp = tmp /. Eps[x_, y_, v_, z_] FV[s1_, v_] FV[s2_, z_] ->
    Eps[x, y, v, z] If[s2 >= s1,  FV[s1, v] FV[s2, z], - FV[s1, z] FV[s2, v] ];

tmp = tmp /. SP[x_,y_] -> sprod[x,y];   
tmp = ToLaTeXForm[tmp];  

tmp = StringReplace[tmp,"Generation$\$$"->"g_"]; 
tmp = StringReplace[tmp,"Gluon$\$$"->"b_"];
If[DiagType == 2,
   tmp = StringReplace[tmp,"mu$\$$1"->"\\mu"];
   tmp = StringReplace[tmp,"mu$\$$2"->"\\nu"]
];
tmp = StringReplace[tmp, " ^" -> "^"]; 

tmp = StringReplace[tmp,"{}"->""];
tmp = StringReplace[tmp, "\\text{" ~~ Shortest[x__] ~~ "}" -> x];
tmp = StringReplace[tmp, "\\text" -> ""];

tmp = StringReplace[tmp,"hlambda"->"\\lambda"];
tmp = StringReplace[tmp,"{vev}"->"{v}"];
tmp = StringReplace[tmp,"vev"->"v"];
tmp = StringReplace[tmp,"MH"->"M_h"];
tmp = StringReplace[tmp,"MW"->"M_W"];
tmp = StringReplace[tmp,"MZ"->"M_Z"];
tmp = StringReplace[tmp,"G_s"->"{\\bar g}_s"];
(*
tmp = StringReplace[tmp,"g1norm"->"Z_{g^{\\prime}}"];
tmp = StringReplace[tmp,"gwnorm"->"Z_{g}"];
tmp = StringReplace[tmp,"gsnorm"->"Z_{g_s}"];
tmp = StringReplace[tmp,"G0norm"->"Z_{G^0}"];
tmp = StringReplace[tmp,"GPnorm"->"Z_{G^+}"];
tmp = StringReplace[tmp,"Hnorm"->"Z_h"];
tmp = StringReplace[tmp,"Wnorm"->"Z_W"];
tmp = StringReplace[tmp,"Gnorm"->"Z_{gl}"];
tmp = StringReplace[tmp,"AZnorm11"->"Z_{\\gamma Z}^{11}"];
tmp = StringReplace[tmp,"AZnorm12"->"Z_{\\gamma Z}^{12}"];
tmp = StringReplace[tmp,"AZnorm21"->"Z_{\\gamma Z}^{21}"];
tmp = StringReplace[tmp,"AZnorm22"->"Z_{\\gamma Z}^{22}"];
*)   
tmp = StringReplace[tmp,"Z_{g'}"->"Z_{g^{\\prime}}"];
tmp = StringReplace[tmp,"Z_{gs}"->"Z_{g_s}"];
tmp = StringReplace[tmp,"Z_{G0}"->"Z_{G^0}"];
tmp = StringReplace[tmp,"Z_{G+}"->"Z_{G^+}"];

tmp = StringReplace[tmp,"Z_{$\\gamma $Z}_{1,1}" -> "Z_{\\gamma Z}^{11}"];
tmp = StringReplace[tmp,"Z_{$\\gamma $Z}_{1,2}" -> "Z_{\\gamma Z}^{12}"];
tmp = StringReplace[tmp,"Z_{$\\gamma $Z}_{2,1}" -> "Z_{\\gamma Z}^{21}"];
tmp = StringReplace[tmp,"Z_{$\\gamma $Z}_{2,2}" -> "Z_{\\gamma Z}^{22}"];

tmp = StringReplace[tmp," ^"->"^"];
tmp = StringReplace[tmp,"^ "->"^"];
tmp = StringReplace[tmp," _"->"_"];
tmp = StringReplace[tmp,"_ "->"_"];

tmp = StringReplace[tmp, "AZnorm(" ~~ Shortest[x__] ~~ "," ~~ Shortest[y__] ~~ ")"->
                         "Z_{\\gamma Z}^{" ~~ x ~~ y ~~ "}"];
tmp = StringReplace[tmp, "AZnorm_{" ~~ Shortest[x__] ~~ "," ~~ Shortest[y__] ~~ "}"->
                         "Z_{\\gamma Z}^{" ~~ x ~~ y ~~ "}"];

tmp = StringReplace[tmp, "Z_{\\gamma Z}^{11}^"-> "\\left(Z_{\\gamma Z}^{11}\\right)^"];
tmp = StringReplace[tmp, "Z_{\\gamma Z}^{12}^"-> "\\left(Z_{\\gamma Z}^{12}\\right)^"];
tmp = StringReplace[tmp, "Z_{\\gamma Z}^{21}^"-> "\\left(Z_{\\gamma Z}^{21}\\right)^"];
tmp = StringReplace[tmp, "Z_{\\gamma Z}^{22}^"-> "\\left(Z_{\\gamma Z}^{22}\\right)^"];

tmp = StringReplace[tmp,"AAAB"->"{\\bar g}{}"];
tmp = StringReplace[tmp,"AAAA^2"->"{\\bar g}^{\\prime 2}"];
tmp = StringReplace[tmp,"AAAA^3"->"{\\bar g}^{\\prime 3}"];
tmp = StringReplace[tmp,"AAAA^4"->"{\\bar g}^{\\prime 4}"];
tmp = StringReplace[tmp,"AAAA"->"{\\bar g}^\\prime"];

tmp = StringReplace[tmp,"{Lambda}"->" "];
(* tmp = StringReplace[tmp,"{Lambda}"->"\\frac{1}{\\Lambda^2}"]; *)

tmp = StringReplace[tmp,"  "->" "];
tmp = StringReplace[tmp," _"->"_"];
tmp = StringReplace[tmp,"_ "->"_"];

tmp = StringReplace[tmp,"P_-" -> "P_L"];
tmp = StringReplace[tmp,"P_+" -> "P_R"];

tmp = StringReplace[tmp,"+"->" + "];
tmp = StringReplace[tmp,"-"->" - "];

tmp = StringReplace[tmp, "alpha " ~~ x_ -> "alpha_" ~~ x];
tmp = StringReplace[tmp, "beta " ~~ x_ -> "beta_" ~~ x];
tmp = StringReplace[tmp, "gamma " ~~ x_ -> "gamma_" ~~ x];

tmp = StringReplace[tmp, "fml\\left(" ~~ Shortest[x__] ~~ "," ~~ Shortest[y__] ~~ "\\right)"->
		    "m_{l_{" ~~ x ~~ "}} \\delta_{" ~~ x ~~ " " ~~ y ~~ "}" ];
tmp = StringReplace[tmp, "fmv\\left(" ~~ Shortest[x__] ~~ "," ~~ Shortest[y__] ~~ "\\right)"->
		    "m_{\\nu_{" ~~ x ~~ "}} \\delta_{" ~~ x ~~ " " ~~ y ~~ "}" ];
tmp = StringReplace[tmp, "fmu\\left(" ~~ Shortest[x__] ~~ "," ~~ Shortest[y__] ~~ "\\right)"->
		    "m_{u_{" ~~ x ~~ "}} \\delta_{" ~~ x ~~ " " ~~ y ~~ "}" ];
tmp = StringReplace[tmp, "fmd\\left(" ~~ Shortest[x__] ~~ "," ~~ Shortest[y__] ~~ "\\right)"->
		    "m_{d_{" ~~ x ~~ "}} \\delta_{" ~~ x ~~ " " ~~ y ~~ "}" ];

tmp = StringReplace[tmp, "{\\nu }" -> "{\\nu}"];
tmp = StringReplace[tmp, "m^{\\nu}_{" ~~ Shortest[x__] ~~ "," ~~ Shortest[y__] ~~"}" ->
		    "m_{\\nu_{" ~~ x ~~ "}} \\delta_{" ~~ x ~~ " " ~~ y ~~ "}" ];
tmp = StringReplace[tmp, "m^l_{" ~~ Shortest[x__] ~~ "," ~~ Shortest[y__] ~~"}" ->
		    "m_{l_{" ~~ x ~~ "}} \\delta_{" ~~ x ~~ " " ~~ y ~~ "}" ];
tmp = StringReplace[tmp, "m^u_{" ~~ Shortest[x__] ~~ "," ~~ Shortest[y__] ~~"}" ->
		    "m_{u_{" ~~ x ~~ "}} \\delta_{" ~~ x ~~ " " ~~ y ~~ "}" ];
tmp = StringReplace[tmp, "m^d_{" ~~ Shortest[x__] ~~ "," ~~ Shortest[y__] ~~"}" ->
		    "m_{d_{" ~~ x ~~ "}} \\delta_{" ~~ x ~~ " " ~~ y ~~ "}" ];

tmp = StringReplace[tmp, "{Subsuperscript}[p," ~~ x_ ~~ "," ~~ Shortest[y__] ~~ "]" ->
                         "p_{" ~~ x ~~ "}^{" ~~ y ~~ "}"];
tmp = StringReplace[tmp, "{Subsuperscript}\\left[p," ~~ x_ ~~ "," ~~ Shortest[y__] ~~ "\\right]" ->
                         "p_{" ~~ x ~~ "}^{" ~~ y ~~ "}"];
tmp = StringReplace[tmp, "sprod$\\$$" ~~ Shortest[z___] ~~ "(" ~~ Shortest[x__] ~~ "," ~~ Shortest[y__] ~~ ")" ->
                         "p_{" ~~ x ~~ "}\\nobreak\\cdot\\nobreak{}p_{" ~~ y ~~ "}"];

tmp = StringReplace[tmp,"$"->" "]; 

tmp = StringReplace[tmp, "m_" ~~ x_ ~~ "_{" ~~ Shortest[y__] ~~ "}" ->
                         "m_{" ~~ x ~~ "_{" ~~ y ~~ "}}"];

tmp = StringReplace[tmp, "SlashedP(" ~~ Shortest[x__] ~~ ")" ->
                         "\\slashed{p}_{" ~~ x ~~ "}"];
(* tmp = StringReplace[tmp,"SIG" ~~ Shortest[x__] ~~ "\\left(" -> "SIG\\left("];*)
If[DiagType == 2, 
(*   tmp = StringReplace[tmp, "\\gamma^{\\mu}.P_" ~~ Shortest[x_] ~~ "_{s_" ~~ Shortest[s1_] ~~
		       ",s_" ~~ Shortest[s2_] ~~ "} \\gamma^{\\mu}.P_" ~~ Shortest[y_] ~~
		       "_{s_" ~~ Shortest[s3_] ~~ ",s_" ~~ Shortest[s4_] ~~ "}" ->
		       "(\\gamma^{\\mu}.P_" ~~ x ~~ ")_{s_" ~~ s1 ~~ " s_" ~~ s2 ~~
		       "} (\\gamma_{\\mu}.P_" ~~ y ~~ ")_{s_" ~~ s3 ~~ " s_" ~~ s4 ~~ "}"];*)
(* with the replacement in line below all Lorentz indices are always up... *)
   tmp = StringReplace[tmp, "\\gamma^{\\mu}.P_" ~~ Shortest[x_] ~~ "_{s_" ~~ Shortest[s1_] ~~
		       ",s_" ~~ Shortest[s2_] ~~ "}" ->
		       "(\\gamma^{\\mu}.P_" ~~ x ~~ ")_{s_" ~~ s1 ~~ " s_" ~~ s2 ~~ "}"];  
   tmp = StringReplace[tmp, "SIG\\left(\\mu,\\nu,P_" ~~ Shortest[x_] ~~ ",s_" ~~ Shortest[s1_] ~~
		       ",s_" ~~ Shortest[s2_] ~~ "\\right) SIG\\left(\\mu,\\nu,P_" ~~
		       Shortest[y_] ~~ ",s_" ~~ Shortest[s3_] ~~ ",s_" ~~ Shortest[s4_] ~~ "\\right)" ->
		       "(\\sigma^{\\mu\\nu} P_" ~~ x ~~ ")_{s_" ~~ s1 ~~ " s_" ~~ s2 ~~
		       "} (\\sigma_{\\mu\\nu} P_" ~~ y ~~ ")_{s_" ~~ s3 ~~ " s_" ~~ s4 ~~ "}" ]; 
];

tmp = StringReplace[tmp, "{Superscript}[C," ~~ Shortest[x__] ~~ "]" ~~
                                               Shortest[y__] ~~ "^*" ->
                                       "C^{" ~~ x ~~  "*}" ~~ y];
tmp = StringReplace[tmp, "SIG\\left(" ~~ Shortest[x__] ~~ "," ~~ Shortest[y__] ~~ "," ~~
                                         Shortest[z__] ~~ "," ~~ Shortest[v__] ~~ "\\right)" ->
                         "\\sigma^{" ~~ x ~~ " " ~~ y ~~ "} " ~~ z ~~ " "];

tmp = StringReplace[tmp,"Gtilde" -> "{\\widetilde G}"];
tmp = StringReplace[tmp,"Wtilde" -> "{\\widetilde W}"];
tmp = StringReplace[tmp,"Btilde" -> "{\\widetilde B}"];

If [DiagType == 0,
  tmp = StringReplace[tmp,"\\delta_{m_1,m_2}" -> ""];
];
If [DiagType == 0 || DiagType == 1,
  tmp = StringReplace[tmp,"\\delta_{s_1,s_2}" -> ""];
  tmp = StringReplace[tmp,"_{s_1,s_2}" -> ""];
];

tmp = StringReplace[tmp, "P_"~~ Shortest[x_] ~~ "_{s_" ~~ Shortest[s1_] ~~ ",s_" ~~ Shortest[s2_] ~~ "}" ->
                        "(P_"~~ x ~~ ")_{s_" ~~ s1 ~~ " s_" ~~ s2 ~~ "}"];

tmp = StringReplace[tmp, "^{" ~~ Shortest[x__] ~~ "}_{" ~~ Shortest[y__] ~~ "}^*" -> "^{" ~~ x ~~ "*}_{" ~~ y ~~ "}"];

tmp = StringReplace[tmp,"Box" -> "\\Box"];
tmp = StringReplace[tmp,"\\phi" -> "\\varphi"];

tmp = StringReplace[tmp,"T_{" -> "{\\cal T}_{"];

tmp = StringReplace[tmp,"."->" "];
tmp = StringReplace[tmp,","->" "];

If[verbose,
   Print[tmp];
   Print[ ];
];

tmp <> "\n"

]]
(* end of ToLaTex *)
]




SelectFieldLine = Function[{n,x1,y1,x2,y2},
Module[{s,f,line,xs,ys,s1,c1,x3,x4,y3,y4,a1,a2,l},

s = StringJoin["(",ToString[x1],",",ToString[y1],")"];
f = StringJoin["(",ToString[x2],",",ToString[y2],")"];

xs = (x1 + x2)/2;
ys = (y1 + y2)/2;

s1 =   (x2 - x1 - y2 + y1)/Sqrt[2]/Sqrt[(x2-x1)^2 + (y2-y1)^2];
c1 = - (x2 - x1 + y2 - y1)/Sqrt[2]/Sqrt[(x2-x1)^2 + (y2-y1)^2];

x3 = N[xs + 8 c1, {Infinity, 1}];
y3 = N[ys + 8 s1, {Infinity, 1}];
x4 = N[xs - 8 s1, {Infinity, 1}];
y4 = N[ys + 8 c1, {Infinity, 1}];
xs = N[xs, {Infinity, 1}];
ys = N[ys, {Infinity, 1}];

If[Abs[x3]<0.5, x3=0];
If[Abs[y3]<0.5, y3=0];
If[Abs[x4]<0.5, x4=0];
If[Abs[y4]<0.5, y4=0];
If[Abs[xs]<0.5, xs=0];
If[Abs[ys]<0.5, ys=0];

a1 = StringJoin["(",ToString[xs],",",ToString[ys],")(",ToString[x3],",",ToString[y3],")"];
a2 = StringJoin["(",ToString[xs],",",ToString[ys],")(",ToString[x4],",",ToString[y4],")"];

Switch[n,
  1, StringJoin["\\DashLine",s,f,"{3}\n"],
  2, StringJoin["\\DashArrowLine",s,f,"{3}\n"],
  3, StringJoin["\\Gluon",s,f,"{3}{4}\n"],
  4, StringJoin["\\Photon",s,f,"{3}{4}\n"],
  5, StringJoin["\\ArrowLine",f,s,"\n"],
  6, StringJoin["\\ArrowLine",s,f,"\n"],
  7, StringJoin["\\Photon",s,f,"{3}{4}\n\\Line",a1,"\n\\Line",a2,"\n"],
  8, StringJoin["\\ZigZag",s,f,"{3}{4}\n"],
  9, StringJoin["\\ZigZag",s,f,"{3}{4}\n"]
]

]
(* end of SelectFieldLine *)
]



DrawDiagram = Function[{cfile,np,feynrule,DiagType,verbose,vname},
Module[{i,x1,x2,y1,y2,tp,pos,vpos,ysize,yshift,vfile,frule,f,fbar},

frule = feynrule;
       
f = {(*vl,*)l,uq,dq};
fbar =  {(*vlbar,*)lbar,uqbar,dqbar};
       
(* Shall we replace clashing arrows with ^C fields with inverted fermion flow? *)
If [ MemberQ[ f, frule[[1,1,1]] ], 
  frule[[1,1,1]] = ToExpression[ ToString[ frule[[1,1,1]] ] <> "cbar"]; 
];
If [ MemberQ[ fbar, frule[[1,2,1]] ], 
  frule[[1,2,1]] = ToExpression[ StringDrop[ ToString[ frule[[1,2,1]] ], -3] <> "c"]; 
];
If [ MemberQ[ f, frule[[1,3,1]] ], 
  frule[[1,3,1]] = ToExpression[ ToString[ frule[[1,3,1]] ] <> "cbar"]; 
];
If [ np > 3 && MemberQ[ fbar, frule[[1,4,1]] ], 
  frule[[1,4,1]] = ToExpression[ StringDrop[ ToString[ frule[[1,4,1]] ], -3] <> "c"]; 
];
       
x1 = {{5, 40, 75},
      {5, 40, 75, 40},
      {5, 20, 60, 75, 40},
      {5, 20, 60, 75, 20, 60}};

y1 = {{0,  50,  0},
      {50, 100, 50, 0},
      {50, 95,  95, 50, 0},
      {50, 95,  95, 50, 5, 5}};

x2 = {{40, 40, 40},
      {40, 40, 40, 40},
      {40, 40, 40, 40, 40},
      {40, 40, 40, 40, 40, 40}};

y2 = {{0, 0, 0},
      {50, 50, 50, 50},
      {50, 50, 50, 50, 50},
      {50, 50, 50, 50, 50, 50}};

tp = {{"(0,0)[r]", "(45,45)[l]","(80,0)[l]"},
      {"(0,50)[r]","(45,95)[l]","(80,50)[l]","(45,5)[l]"},
      {"(0,50)[r]","(15,90)[r]","(65,90)[l]","(80,50)[l]","(45,5)[l]"},
      {"(0,50)[r]","(15,90)[r]","(65,90)[l]","(80,50)[l]","(15,10)[r]","(65,10)[l]"}};

ysize = {"55)","100)","100)","100)"};
yshift = {"(0,15)\n","(0,50)\n","(0,50)\n","(0,50)\n"};
vpos =  {0,50,50,50};

WriteString[cfile, "\\noindent \\begin{tabular}{cp{1mm}l}\n"];
WriteString[cfile, "\\begin{picture}(90,"~~ysize[[np-2]]~~yshift[[np-2]] ];

For[i=1,i<np+1,i++, 
    pos = Position[ FieldType, frule[[1]][[i]][[1]] ][[1]][[1]]; 
    WriteString[cfile, SelectFieldLine[LineType[[pos]],x1[[np-2,i]],y1[[np-2,i]],x2[[np-2,i]],y2[[np-2,i]]] ];
    WriteString[cfile, "\\Text"~~tp[[np-2,i]]~~"{"~~FieldNameAndIndex[pos,i,DiagType]~~"}\n"];
];

WriteString[cfile, "\\Vertex(40,"~~ToString[ vpos[[np-2]] ]~~"){2}\n"];
WriteString[cfile, "\\end{picture}\n"];
WriteString[cfile, "&&\n"];
WriteString[cfile, "\\begin{minipage}[c]{0.8\\linewidth}\n"];
WriteString[cfile, "\\input vertices/"~~vname[[1]]~~"/"~~vname[[2]]~~".tex\n" ];
WriteString[cfile, "\\end{minipage}\n"];
WriteString[cfile, "\\end{tabular}\n\n"];
WriteString[cfile, "\\bigskip\n\n"];

vfile = OpenWrite[ FileNameJoin[{SMEFT$Path, "output", "latex",
"vertices", vname[[1]], vname[[2]]<>".tex"}] ];

WriteString[vfile, "%\n\\begin{dmath*}\n%\n"];
WriteString[vfile, ToLaTeX[frule[[2]], DiagType, verbose ] ];
WriteString[vfile, "%\n\\end{dmath*}\n%\n"];
Close[vfile];


];
(* end of DrawDiagram *)
];


DrawSector = Function[{cfile,frset,verbose,DiagType,DiagClass},
(* draws diagrams for chosen sector of Lagrangian *)		
Module[{i,j,vlist,vtype,vcount,vname},

vlist = {};
vcount = 0;

For[j=3,j<7,j++,
  For[i=1,i < Length[frset]+1,i++, 
    vtype = Dim6VertexName[ frset[[i,1]] ];
    If[ ! MemberQ[vlist, vtype] && Length[frset[[i]][[1]]] == j && frset[[i]][[2]] =!= 0,
    If[ verbose, Print[vtype] ];
      vname = FileNameJoin[{SMEFT$Path, "output", "latex", "vertices", DiagClass}];
      If [ ! DirectoryQ[vname],  CreateDirectory[vname] ];
      vname = {DiagClass, vtype};
      WriteString[cfile, "%Vertex " <> vtype <> "\n"];
      DrawDiagram[cfile,j,frset[[i]], DiagType, verbose, vname ];
      vlist = AppendTo[vlist, vtype];
      vcount = vcount + 1;
    ];
  ];
];

vcount

]
(* end of DrawSector *)
]




SMEFTToLatex[ OptionsPattern[{ FullDocument -> True, ScreenOutput -> False, Expansion->"none" }] ] :=
(* master routine of latex generation code *)
Module[{i, j, k, frules, cfile, VertCount, vlist, nlist, olist, ghostname, tmp},

tmp = Complement[{ OptionValue[FullDocument], OptionValue[ScreenOutput]}, {True,False}];
If[ tmp =!= {},
  Print["Options FullDocument and ScreenOutput of SMEFTToLatex can be only True or False, value(s) ", tmp, " not allowed, please correct!"];
  Abort[]
];

If [ ! MemberQ[ {"NONE", "SMEFT", "none", "smeft"}, OptionValue[Expansion] ],
  Print["Option Expansion of SMEFTToLatex can be only \"none\" or \"smeft\", please correct!"];
  Abort[]
];
       
vlist = # <> "Vertices" <> If [ OptionValue[Expansion] === "smeft",
"Latex", "" ] & /@ {"LeptonGauge", "LeptonHiggsGauge", "QuarkGauge",
"QuarkHiggsGauge", "QuarkGluon", "GaugeSelf", "GaugeHiggs",
"GluonSelf", "GluonHiggs", "FourLepton", "TwoQuarkTwoLepton",
"FourQuark", "BLViolating"};

If [ OptionValue[Expansion] === "smeft", SMEFTExpandVertices[ ExpOrder-> 1, ExpName-> "Latex"] ];
    
nlist = {"Lepton--gauge", "Lepton--Higgs--gauge", "Quark--gauge",
"Quark--Higgs--gauge", "Quark-gluon", "Gauge self interaction",
"Higgs--gauge", "Gluon self interaction", "Higgs--gluon", "Four
lepton", "Two quark--two lepton", "Four quark", "Baryon and lepton
number violating four fermion"};

olist = {0,0,0,0,1,0,0,1,1,2,2,2,2};

Print[ Style["Generating Latex file with Feynman rules...", Bold] ];
VertCount = 0;

(* open file and create header *)
cfile = OpenWrite[ FileNameJoin[{SMEFT$Path, "output", "latex", "smeft_feynman_rules.tex"}] ];

Print["Latex form of Feynman rules will be stored in file ",
  Style[ FileNameJoin[{SMEFT$Path, "output", "latex", "smeft_feynman_rules.tex"}], Bold]];
Print["Compilation requires breqn.sty and axodraw.sty styles (included)."];
Print["Automatic line breaking of math formulae may not be perfect!"];

WriteString[cfile, "%Automatically generated Feynman rules file for SMEFT\n\n"];
If[OptionValue[FullDocument],
   WriteString[cfile, "\\documentclass[11pt]{article}\n"];
   WriteString[cfile, "\\usepackage{axodraw}\n"];
   WriteString[cfile, "\\usepackage{amsmath, amssymb}\n"];
   WriteString[cfile, "\\usepackage{slashed}\n"];
   WriteString[cfile, "\\usepackage{breqn}\n"];
   WriteString[cfile, "\\setkeys{breqn}{breakdepth={3}}\n"];

   WriteString[cfile, "\\textwidth = 16cm\n\n"];
   WriteString[cfile, "\\textheight = 24cm\n\n"];
   WriteString[cfile, "\\topmargin=-2cm\n"];
   WriteString[cfile, "\\oddsidemargin=0cm\n"];
   WriteString[cfile, "\\evensidemargin=\\oddsidemargin\n"];

(*   WriteString[cfile, "\\def\\text#1{#1}\n\n"];*)
   WriteString[cfile, "\\begin{document}\n\n"];

   WriteString[cfile, "\\appendix\n\n"];

   WriteString[cfile, "\\bigskip\n\\bigskip\n\\bigskip\n"];
   WriteString[cfile, "\\section{SMEFT interaction vertices}\n\n"];

   WriteString[cfile, "{\\bf CAUTION:} interaction vertices printed below are shown including only terms up to mass dimension-6. Interactions proportional to products of dimension-6 Wilson coefficients, even if calculated and included in other output formats (Mathematica, Feynarts, UFO etc.), are too complicated for printout and for manual calculations. If necessary, they can be inspected visually displaying relevant variables (for their list see {\\tt SmeftFR} manual) in the Mathematica notebook.\n"];
];

For[i=1, i < Length[vlist] + 1, i++, 
  WriteString[cfile, "\\bigskip\n\\bigskip\n"];
  WriteString[cfile, "\\subsection{" <> nlist[[i]] <> " vertices}\n\n"];
  frules = Expand[ ToExpression[ vlist[[i]] ] ] /. Lam^k_ -> If[ k > 1, 0, Lam^k ];
(* skip gluon vertices with more than 4 legs *)
  If[ MemberQ[ {"GluonSelfVertices", "GluonSelfVerticesLatex"}, vlist[[i]] ],  
    For[j=1, j < Length[frules] + 1, j++,
      If[ Length[ frules[[j,1]] ] > 4 , frules[[j,2]] = LongExpressionNotDisplayed ];
    ];
  ];
  VertCount = VertCount + DrawSector[cfile,frules,OptionValue[ScreenOutput],olist[[i]],vlist[[i]]]; 
];

If[SMEFT$RxiGaugeStatus,
   ghostname = If [OptionValue[Expansion] === "smeft", "GhostVerticesLatex", "GhostVertices"];
   WriteString[cfile, "\\bigskip\n\\bigskip\n\\bigskip\n"];
   WriteString[cfile, "\\subsection{Ghost vertices}\n\n"];
   frules = Expand[ ToExpression[ghostname] ] /. Lam^k_ -> If[ k > 1, 0, Lam^k ];
   VertCount = VertCount + DrawSector[cfile,frules,OptionValue[ScreenOutput],0,ghostname];
];
    
If[OptionValue[FullDocument],
   WriteString[cfile, "\\end{document}\n"];
];

Close[ FileNameJoin[{SMEFT$Path, "output", "latex", "smeft_feynman_rules.tex"}] ];

Print["Latex output ready, total number of vertices drawn = ", VertCount];

]
(* end of SMEFTLatexGenerator *)
