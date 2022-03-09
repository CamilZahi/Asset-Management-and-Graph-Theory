data cac1;
 set cac;
 array NumVar _numeric_;
 do over NumVar;
 if NumVar=. then NumVar=0;
 end;
run;

proc iml;
use cac1;
read all var _NUM_ into adjmat[colname=actifs];
close cac1;
use price_no_date;
read all var _NUM_ into price [colname=actifs];
close price_no_date;
adj=j(ncol(adjmat),ncol(adjmat),0);
start vec(mat);
free mat1;
do z=1 to ncol(mat);
 mat1=mat1//mat[,z];
 end;
return(mat1);
finish vec;
start vechh(mat);
free mat2;
do z=1 to nrow(mat);
 mat2=mat2||mat[z,];
 end;
return(mat2);
finish vechh;
prix_vente_total=1000;/*initialisation de notre premi?re somme 
d'investissement. Elle sera r??valu?e ? chaque p?riode*/
prix_vente_totalCAC=1000;
free valorisa;
Max=floor((nrow(adjmat)-132)/44);
valorisa=j(1,4,1);
Do h=0 to Max;
Do f=1 to ncol(adjmat);
 Do g=f to ncol(adjmat);
m=1;
a1=adjmat[1+2*22*h:132+2*22*h,f];
a2=adjmat[1+2*22*h:132+2*22*h,g];
d1=ncol(a1);
d2=ncol(a2);
a=a1||a2;
c=a`*a/nrow(a);
c1=a1`*a1/nrow(a1);
c2=a2`*a2/nrow(a2);
c1=c[1:d1,1:d1];
c2=c[d1+1:nrow(c),d1+1:ncol(c)];
r1=diag(vecdiag(c1)##-0.5);
r2=diag(vecdiag(c2)##-0.5);
r11=(r1)*c1*(r1);
r22=(r2)*c2*(r2);
QhP=j(m+1,3,0);
QhM=j(m+1,3,0);
Icont=j(m+1,d1*2+1,0);/*Individual contribution negative lags*/
Icontm=j(m+1,d2*2+1,0);/*Individual contribution positive lags*/
free qhle qhla vla1 vle1 hla1 hle1 Icor1 Icor2;
compt=1;
do k=0 to m;
k1=m-k;
/*leads of a2, k<0 : a1 drives a2*/
 leads=r1*(a1[1:nrow(a)-k1,]`*a2[1+k1:nrow(a),]/nrow(a1))*r2;
 Vle=vec(leads);
 var=sqrt ((nrow(a)-abs(k1)) / (nrow(a)*nrow(a))) ;
 *Icor1=Icor1//(-k1||vechh(leads)||-2*var||2*var);
 *Vle1=Vle1//Vle;
 QhM[compt,1]=-k1;
 QhM[compt,2]=nrow(a)*((vle`*inv(r22@r11)*vle)*nrow(a)/(nrow(a)-
abs(k1)));
 QhM[compt,3]=1-probchi(abs(QhM[compt,2]),d1*d2);
M1=j(nrow(vle),nrow(vle),0);
 M1[1,1]=1;
 vlec1=M1*vle;
Icont[compt,1]=-k1;
 
Icont[compt,2]=nrow(a)*((vlec1`*inv(r22@r11)*vlec1)*nrow(a)/(nrow(a)-
abs(k1)));
 Icont[compt,2+1]=1-
probchi(abs(Icont[compt,2]),d2);
/*lags of a2, k>0 : a2 drives a1*/
 lags=r1*(a1[1+k:nrow(a),]`*a2[1:nrow(a)-k,]/nrow(a))*r2;
 vla=vec(lags);
 var=sqrt ((nrow(a)-abs(k)) / (nrow(a)*nrow(a))) ;
 Icor2=Icor2//(k||vechh(lags)||-2*var||2*var);
 QhP[compt,1]=k;
 QhP[compt,2]=nrow(a)*((vla`*inv(r22@r11)*vla)*nrow(a)/(nrow(a)-
abs(k)));
 QhP[compt,3]=1-probchi(abs(QhP[compt,2]),d1*d2);
 M2=j(nrow(vla),nrow(vla),0);
M2[1,1]=1;
vlags1=M2*vla;
 Icontm[compt,1]=k;
Icontm[compt,2]=nrow(a)*((vlags1`*inv(r22@r11)*vlags1)*nrow(a)/(nrow(a)-
abs(k)));
 Icontm[compt,2+1]=1-
probchi(abs(Icontm[compt,2]),d1);
compt=compt+1;
end;
plus='k';
plus=plus//concat('X1',char(1,1.0))//'Pvalue';
plusm='k';
plusm=plusm//concat('X2',char(1,1.0))//'Pvalue';
/*Individual correlations*/
icor=icor1//icor2[2:nrow(icor2),];
col1=j(1,1,'rho');
col2=j(1,1,'X1');
col3=j(1,1,'X2');
*free col4;
*col4=col4//concat('rho',col2[1],char(1,1),'_',col3[1],char(1,1));
*col4='k'//col4//'Down95'//'Up95';
*print cal4;
Icause=j(1,ncol(Icont)-1,0);
Icause[1]=sum(Icont[1:nrow(Icont)-1,2]);
Icause[2]=1-probchi(abs(Icause[1]),m*d2);
plus1=plus[2:nrow(plus)];
plusm1=plusm[2:nrow(plusm)];
Icausem=j(1,ncol(Icontm)-1,0);
Icausem[1]=sum(Icontm[2:nrow(Icont),2]);
Icausem[2]=1-probchi(abs(Icausem[1]),m*d1);
*******************;
lagmax=m;
pccmX1toX2=j(1,4,0);
pccmX2toX1=j(1,4,0);
do lag=1 to lagmax;
free Xa;
/*J model*/
do i=1 to lag;
S1=I(nrow(a1)-i);
S2=j(nrow(a1),nrow(a1)-nrow(s1),0);
s3=j(nrow(s2)-nrow(s1),nrow(s1),0);
S=s3//s1||s2;
Xa=Xa||S*a1;
end;
Xab=Xa[lag+1:nrow(Xa),];
Yab=a2[lag+1:nrow(Xa),];
beta=inv(Xab`*Xab)*Xab`*Yab;
r=Yab-Xab*beta;
sygmaJ=det((r`*r)/(nrow(a2)));
free Xa;
/*J-1 model*/
do i=1 to lag-1;
S1=I(nrow(a1)-i);
S2=j(nrow(a1),nrow(a1)-nrow(s1),0);
s3=j(nrow(s2)-nrow(s1),nrow(s1),0);
S=s3//s1||s2;
Xa=Xa||S*a1;
end;
if lag=1 then do;
Yab=a2[lag+1:nrow(a2),];
sygmaJ1=det((Yab`*Yab)/(nrow(a2)));
end;
else do;
Xab=Xa[lag+1:nrow(Xa),];
Yab=a2[lag+1:nrow(Xa),];
beta1=inv(Xab`*Xab)*Xab`*Yab;
r=Yab-Xab*beta1;
sygmaJ1=det((r`*r)/(nrow(a2)));
end;
test=-(nrow(a1)-lag*d1*d2-1-0.5)*log(SygmaJ/SygmaJ1);
if test>=0 then test=test;else test=0;
ptest=1-probchi(abs(test),d1*d2);
pccmX1toX2[lag,1]=-lag;
pccmX1toX2[lag,2]=test;
pccmX1toX2[lag,3]=ptest;
end;
mattrib pccmX1toX2[Colname={'k','LR test','P-value'} label='Test based 
on Partial Cross-Correlation Function at k'];
call sort(pccmX1toX2,1);
/*********From X2 to X1********/
do lag=1 to lagmax;
free Xa;
/*J model*/
do i=1 to lag;
S1=I(nrow(a2)-i);
S2=j(nrow(a2),nrow(a2)-nrow(s1),0);
s3=j(nrow(s2)-nrow(s1),nrow(s1),0);
S=s3//s1||s2;
Xa=Xa||S*a2;
end;
Xab=Xa[lag+1:nrow(Xa),];
Yab=a1[lag+1:nrow(Xa),];
beta=inv(Xab`*Xab)*Xab`*Yab;
r=Yab-Xab*beta;
sygmaJ=det((r`*r)/(nrow(a2)));
free Xa;
/*J-1 model*/
do i=1 to lag-1;
S1=I(nrow(a2)-i);
S2=j(nrow(a2),nrow(a2)-nrow(s1),0);
s3=j(nrow(s2)-nrow(s1),nrow(s1),0);
S=s3//s1||s2;
Xa=Xa||S*a2;
end;
if lag=1 then do;
Yab=a1[lag+1:nrow(a1),];
sygmaJ1=det((Yab`*Yab)/(nrow(a1)));
end;
else do;
Xab=Xa[lag+1:nrow(Xa),];
Yab=a1[lag+1:nrow(Xa),];
beta1=inv(Xab`*Xab)*Xab`*Yab;
r=Yab-Xab*beta1;
sygmaJ1=det((r`*r)/(nrow(a2)));
end;
*print lag sygmaJ sygmaJ1 (log(sygmaJ/sygmaJ1));
test1=-(nrow(a1)-lag*d1*d2-1-0.5)*log(SygmaJ/SygmaJ1);
if test1>=0 then test=test;else test=0;
ptest=1-probchi(abs(test1),d1*d2);
pccmX2toX1[lag,1]=lag;
pccmX2toX1[lag,2]=test1;
pccmX2toX1[lag,3]=ptest;
end;
pccm=pccmX1toX2//{0 . . .}//pccmX2toX1;
mattrib pccm[Colname={'k','LR test','P-value', 'Sign.(5%)'} label='Test based 
on Partial Cross-Correlation Function at lag k'];
do i=1 to nrow(pccm);
if pccm[i,3]<=0.05 & pccm[i,3]^=. then pccm[i,4]=11111;
else pccm[i,4]=.;
end;
/*Causality Using an LR test: X1 on X2*/
free Xa;
/*J model*/
do i=1 to lagmax;
S1=I(nrow(a1)-i);
S2=j(nrow(a1),nrow(a1)-nrow(s1),0);
s3=j(nrow(s2)-nrow(s1),nrow(s1),0);
S=s3//s1||s2;
Xa=Xa||S*a1;
end;
Xab=Xa[lagmax+1:nrow(Xa),];
Yab=a2[lagmax+1:nrow(Xa),];
beta=inv(Xab`*Xab)*Xab`*Yab;
r=Yab-Xab*beta;
sygmaJ=det((r`*r)/(nrow(a2)));
free Xa;
/*J-1 model*/
Yab=a2[lagmax+1:nrow(a2),];
sygmaJ1=det((Yab`*Yab)/(nrow(a2)));
LrcausX1toX2=-(nrow(a1)-lagmax*ncol(a1)*ncol(a2)-1-0.5)*log(SygmaJ/SygmaJ1);
PLrcausX1toX2=1-probchi(abs(LrcausX1toX2),lagmax*ncol(a1)*ncol(a2)); 
/*Causality Using an LR test: X2 on X1*/
free Xa;
do i=1 to lagmax;
S1=I(nrow(a2)-i);
S2=j(nrow(a2),nrow(a2)-nrow(s1),0);
s3=j(nrow(s2)-nrow(s1),nrow(s1),0);
S=s3//s1||s2;
Xa=Xa||S*a2;
end;
Xab=Xa[lagmax+1:nrow(Xa),];
Yab=a1[lagmax+1:nrow(Xa),];
beta=inv(Xab`*Xab)*Xab`*Yab;
r=Yab-Xab*beta;
sygmaJ=det((r`*r)/(nrow(a2)));
free Xa;
/*J-1 model*/
Yab=a1[lagmax+1:nrow(a1),];
sygmaJ1=det((Yab`*Yab)/(nrow(a1)));
LrcausX2toX1=-(nrow(a1)-lagmax*ncol(a1)*ncol(a2)-1-0.5)*log(SygmaJ/SygmaJ1);
PLrcausX2toX1=1-probchi(abs(LrcausX2toX1),lagmax*ncol(a1)*ncol(a2));
LRcause=(LrcausX1toX2||PLrcausX1toX2)//(LrcausX2toX1||PLrcausX2toX1);
mattrib Lrcause[rowname={'Non-Causality from X1 to X2:','Non-Causality from X2 
to X1:'} colname={'LR test','P-Value'}
label='LR tests for non-Causality'];
/*****************************************************************************
**************************************/
mattrib icor[colname=col4 label='Individual correlations Corr[X1(t),X2(t-k)]';
mattrib icont[colname=plus label='Individual contribution based on 
Corr[X1i(t),X2(t-k)]'];
mattrib icontm[colname=plusm label='Individual contribution based on 
Corr[X1(t),X2i(t-k)]'];
mattrib icause[colname=plus1 label='Individual contribution to causality on 
X2'];
mattrib icausem[colname=plusm1 label='Individual contribution to causality on 
X1'];
mattrib Qhm [colname={'k','Qhr','Palue','Bin'} label='El Himdi & Roy, based on
Cov[X1(t),X2(t-k)]'];
mattrib QhP [colname={'k','Qhr','Palue','Bin'} label='El Himdi & Roy, based on
Cov[X1(t),X2(t-k)]'];
eLhr=Qhm//QhP[2:nrow(Qhp),];
eLHr1=eLhr||j(nrow(eLhr),1,.);
do i=1 to nrow(eLhr);
if eLhr1[i,3]<=0.05 then eLhr1[i,4]=11111;
else tt=1;
end;
mattrib eLhr1 [colname={'k','Qhr','Pvalue','Sign.5%'} label='Test based on 
Cov[X1(t),X2(t-k)]'];
Qhr=sum(Qhm[,2])+sum(Qhp[2:nrow(Qhp),2]);/*Overall dependance*/
Qhrm=sum(Qhm[1:nrow(Qhm)-1,2]);/*Causality from X1 to X2*/
Qhrp=sum(Qhp[2:nrow(Qhp),2]);/*Causality from X2 to X1*/
tab=(Qhr||1-probchi(abs(Qhr),(2*M+1)*d1*d2))//(Qhrm||1-
probchi(abs(Qhrm),M*d1*d2))//(Qhrp||1-probchi(abs(Qhrp),M*d1*d2));
mattrib tab[colname={'Qhr','Pvalue'} rowname={'Overall Dependence: ','NonCausality from X1 to X2:','Non-Causality from X2 to X1:'} label='NonCorrelation & Non-Causality'];
 if elhr1[1,3]<0.1 then adj[f,g]=1;
 else adj[f,g]=0;
 if elhr1[3,3]<0.1 then adj[g,f]=1;
 else adj[g,f]=0;
end;
end;
*********matrice d'adjacence retouch?e*****;
DO i=1 to ncol(adj);
If adj[i,i]=1 then adj[i,i]=0;
else citron=8;
end;
adjacence=adj[2:28, 2:28];
*********Coefficient IN*****;
un=j(nrow(adjacence),1,1);
din=adjacence`*un;
nom=vecdiag((adjacence`*adjacence*adjacence));
CCii=vecdiag((adjacence`*adjacence*adjacence))/(din#(din-1));
*********gestion du portefeuille*****;
Do y=1 to nrow(CCii);
If CCii[y,1]<= 0.166 then prix = prix || price[,y];
else melon=6;
end;
alloc=divide(prix_vente_total,ncol(prix));/*achat des actifs - dans notre 
ptf*/
quantite_d_achat=divide(alloc,prix[132+2*22*h,]);
allocCAC=divide(prix_vente_totalCAC,ncol(price));/*achat des actifs - avec 
tous les actifs du cac*/
quantite_d_achatCAC=divide(allocCAC,price[132+2*22*h,]);
prix_quon_prendT=prix[132+2*22+2*22*h,]`;/*periode de gestion --> on vend 
notre ptf*/
prix_vente_total= quantite_d_achat*prix_quon_prendT;
prix_quon_prendTCAC=price[132+2*22+2*22*h,]`;
prix_vente_totalCAC= quantite_d_achatCAC*prix_quon_prendTCAC;
Diff=prix_vente_total-prix_vente_totalCAC;
valo=(h||prix_vente_total||prix_vente_totalCAC||diff);
valorisa=(valorisa//valo);
end;
valorisaFINAL=valorisa[2:nrow(valorisa), 1:ncol(valorisa)];
Notre_valorisation = valorisaFINAL[,2]; 
CAC_valorisation = valorisaFINAL[,3];
Notre_valorisation = (1000//Notre_valorisation); 
CAC_valorisation = (1000//CAC_valorisation); 
/*Calcul des rendements par p?riode*/
max_boucle=nrow(CAC_valorisation)-1;
renta_ptf=j(max_boucle,1,0);
renta_CAC=j(max_boucle,1,0);
DO i=1 to max_boucle;
renta_ptf[i,1]=((Notre_valorisation[i+1,1]/Notre_valorisation[i,1])-1);
renta_cac[i,1]=((CAC_valorisation[i+1,1]/CAC_valorisation[i,1])-1);
end;
/*les esp?rances de rendements*/
Rp=(1/nrow(renta_ptf))*sum(renta_ptf);
Mar=(1/nrow(renta_cac))*sum(renta_cac);
/*Sortino*/
index_dsd=loc(renta_ptf<Mar);
rend_inf_Mar= renta_ptf[index_dsd];
Nbr_rend_Mar = nrow(rend_inf_Mar);
dsd=sqrt((1/Nbr_rend_Mar)*(sum((rend_inf_Mar-Mar)##2)));
Sortino = (Rp-Mar)/dsd;
print dsd Sortino;
/*Tracking error*/
el_1=(sum((renta_ptf-renta_cac)##2))/nrow(sum((renta_ptf-renta_cac)##2));
el_2=(sum((renta_ptf-renta_cac))/nrow(renta_ptf))##2;
TrackingError= sqrt(el_1-el_2);
print TrackingError;
