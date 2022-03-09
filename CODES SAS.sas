proc iml;
/*On stock notre table dans une matrice*/
use Cac;
read all var _NUM_ into adjmat[colname=actifs];
close cac;
/*Définition de la matrice d'adjacence*/
adj=j(ncol(adjmat),ncol(adjmat),0);

Do h=0 to 2;
/*On balaye les actifs 2 à 2 */
Do i=1 to ncol(adjmat);
    Do j=i to ncol(adjmat);
*mat=adjmat[,i]||adjmat[,j];
mat=adjmat[1+2*22*h:132+2*22*h,i]||adjmat[1+2*22*h:132+2*22*h,j];
create toto from mat;
    append from mat;
    close toto;
    /* On va calcule les corrélations croisés des actifs 2 à 2 */
    rname={'Corr[X1(t-1),X2(t)]','  ','Corr[X2(t-1),X1(t)]'};
    submit ;
    proc timeseries data=toto outcrosscorr=un;
    crossvar col1 col2; /*Ici, on cherche à mettre les colonnes i et j (qui était avant les 2 actifs*/
    crosscorr ccf ccfprob ccf2std lag / nlag=1;
    run;
    endsubmit;
    /*On construit ensuite la matrice avec les corrélations croisés et les pvalues*/
    use un;
    read all;
    close un;
    resu1=(lag[1:3])||(CCF[1:3])||(CCFPROB[1:3]);
    mattrib resu1[rowname=rname colname={'Lag','CCf','Prob'}
    label='CCF proc time series'];
    /*Remplissage de la matrice d'adjacence*/
    if resu1[1,3]<0.05 then adj[i,j]=1;
    else adj[i,j]=0;
    if resu1[3,3]<0.05 then adj[j,i]=1;
    else adj[j,i]=0;
end;
end;

adjFinal=  adj[2:28, 2:28];
Ds_calculation=j(1,ncol(adjFinal),1);
Ds_in = Ds_calculation * adjFinal;
Numerateur=adjFinal`*adjFinal*adjFinal;
Dsless=j(1,ncol(adjFinal),1);
Ds_inless=Ds_in-Dsless;
Denominateur=Ds_in`*Ds_inless;
Numii = vecdiag(Numerateur);
Denii = vecdiag(denominateur);
/*calcul final de chaque coefficient, ça donne une matrice (27,1)*/
CCii=Numii/Denii;
print CCii;

end;
