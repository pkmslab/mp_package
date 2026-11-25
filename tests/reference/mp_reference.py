# Independent Python reference for predict_mp() (R/predict-mp.R).
# Re-implements the equations in MP_Mechanistic_Models_2.docx and runs the
# bundled XGBoost model on five illustrative drugs. The expected values in
# tests/testthat/test-predict-mp.R were taken from this output.
# Needs: pip install xgboost numpy.
# Note: the ML output printed here is 10^raw, i.e. it assumes a log10 target.
import math, numpy as np, xgboost as xgb
M="/mnt/user-data/attachments/67e18b7e-471a-585f-b00d-9c0675470ec5/e0e094bf-final_MP_model.json"
b=xgb.Booster(); b.load_model("inst/extdata/finlal_MP_model.json")
def run(MW,LogP,PSA,HBD,HBA,pa,pb,Fu,LogD,typ,bd):
    nan=float('nan'); pa=nan if pa is None else pa; pb=nan if pb is None else pb
    f=lambda pH:1/(1+(0 if pa!=pa else 10**(pH-pa))+(0 if pb!=pb else 10**(pb-pH)))
    fp,fm=f(7.4),f(7.0)
    Fum=Fu**.448/(.000694**.448+Fu**.448); mupu=fm/fp
    d70=10**LogD*mupu; l70=-.88+1.29*math.log10(d70); Dm=10**l70
    cs=-3.912-.015*PSA+3.367*math.log10(MW)-.164*(LogP-LogD); cr=2.793+.179*LogP-.132*HBD
    clr=10**cs/10**cr; fmm=1/(.955/Fum+.045*Dm)
    pd=Fu*mupu*(.955/Fum+.045*Dm)
    lpd=math.exp(.03+2.28*math.log(mupu)+.89*math.log(Fu)+.51*math.log(.955/Fum+.045*Dm))
    kos=clr*Fu/fmm
    # meskin
    ml=float('nan')
    if typ!='neutral':
        eq=None;pk=None
        if typ=='acid': eq='acid'
        elif typ=='base': eq='base';pk=pb
        else:
            ac=abs(pa-7)<abs(pb-7)
            if pa>pb: eq,pk=('acid',None) if ac else ('base',pb)
            else: eq,pk=('base',pa) if ac else ('acid',None)
        ml=10**(2.068-.162*MW**.5-.185*LogP) if eq=='acid' else 10**(.265-.153*LogP-.128*(7.4-pk))
    names=b.feature_names
    d=dict(MW=MW,PSA=PSA,HBD=HBD,HBA=HBA,LogP=LogP,Fu=Fu,pKa_acid=pa,pKa_base=pb,MuPu=mupu,LogD74=LogD,
      BDDCS_1=float(bd==1),BDDCS_2=float(bd==2),BDDCS_3=float(bd==3),BDDCS_4=float(bd==4),
      ion_class_neutral=float(typ=='neutral'),ion_class_acid=float(typ=='acid'),ion_class_base=float(typ=='base'),ion_class_ampholyte=float(typ=='ampholyte'))
    raw=b.predict(xgb.DMatrix(np.array([[d[n] for n in names]]),feature_names=names))[0]
    return dict(MuPu=mupu,ML_raw=float(raw),ML_log10=10**raw,PD=pd,LogPD=lpd,Kosh=kos,Meskin=ml)
# test drugs (illustrative round numbers, not literature values)
cases={
 "acid":  (180.16,1.31,63.6,1,3,3.5,None,0.20,1.31+math.log10(1/(1+10**(7.4-3.5))),"acid",3),
 "base":  (285.4,2.9,30.0,1,3,None,9.0,0.60,2.9+math.log10(1/(1+10**(9.0-7.4))),"base",1),
 "neutral":(194.19,-0.07,58.4,0,6,None,None,0.70,-0.07,"neutral",1),
 "amph1": (350.0,2.0,80.0,2,5,8.0,4.0,0.40,1.9,"ampholyte",2),
 "amph2": (350.0,2.0,80.0,2,5,5.0,9.5,0.40,0.5,"ampholyte",2),
}
for k,v in cases.items(): print(k,{a:round(c,5) for a,c in run(*v).items()})
