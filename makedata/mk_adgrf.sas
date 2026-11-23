/******************************************************************************/
/* PROD PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/makedata/production/mk_adgrf.sas
/* WORK PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/makedata/working/mk_adgrf.sas
/* 
/* RISK ASSESSMENT: MODERATE
/* DOCUMENT ID:     PROG.10083.v4
/* 
/* PURPOSE:         Create ADSL for use with grf library in R
/* 
/* SOURCE PRGM:     none
/* INPUT:           derived.m01_adsl, derived.m061_adttmace, derived.m081_adttae_final, 
/*                  derived.m19_adbiomrkrs, derived.m02_riskscore
/* OUTPUT:          derived.adgrf
/*
/* MACROS USED:     none
/* EXEMPTIONS:      none
/* 
/* AUTHOR:          Emma Davies Smith
/* CREATION DATE:   18JUL2025
/* 
/* NOTES:           Competing events recoded to censored
/*                  Using first PRIMARY MACE and CONFIRMED DIABETES
/*                  Converted to adgrf.RDS using sister program mk_adgrf.R
/*     
/* MODIFICATIONS:   none
/******************************************************************************/

/*******************************************************************************
                                    M01_ADSL                                
*******************************************************************************/

data adsl;
    set derived.m01_adsl;
    keep
        /* IDENTIFIERS */
        patid
        rxl
        rxln
        ascvd
        /* DEMOGRAPHICS */
        age
        agegp3
        natalsex
        racegp
        cntryg1
        gbdgp1
        /* CARDIOMETABOLIC */
        cmdiabfl
        cmhypfl
        bmi
        bmigp6
        schol
        diabp
        diabfl
        egfr_eq21
        fhcvd
        shdl
        fgluc
        hypfl
        sldlc
        metS
        prediab
        sysbp
        strig
        wstcir
        ms_wcir
        /* HIV1 */
        acat
        acat2
        abcfl
        hiv1rnaq
        hiv1rnagp2
        nadirgp2
        cd4
        cd4gp3
        arvdur
        arvdurgp
        /* LIFESTYLE */
        cigstgp
        cigncf
        subsgp1
        subsgp2
        ;
    run;

proc contents data=adsl; run;

proc print data=adsl(obs=10);
    title 'Preview selected adsl variables';
run;

/*******************************************************************************
                            DERIVED COVARIATES                                
*******************************************************************************/
data adsl_derived;
    /* ORGANIZE */
    retain  patid rxl rxln pce pcefl age agegp agegp3 natalsex racegp cntryg1 gbdgp1 gbdrace
            cmdiabfl cmhypfl bmi bmigp6 schol diabp diabfl egfr_eq21 egfrfl fhcvd fhcvd2 shdl hdlfl fgluc 
            fglucfl hypfl sldlc ldlfl metS metS2 prediab prediab2 sysbp strig trigfl wstcir ms_wcir
            acat acat2 abcfl instifl nadirgp2 cd4 cd4gp3 arvdur arvdurgp hiv1rnaq hiv1rnagp2
            cigstgp cigncf subsgp1 subsgp2;
    set adsl(rename=(cmdiabfl=_cmdiabfl cmhypfl=_cmhypfl diabfl=_diabfl hypfl=_hypfl ascvd=pce
                     abcfl=_abcfl nadirgp2=_nadirgp2 arvdurgp=_arvdurgp cigstgp=_cigstgp));

    /* IDENTIFIERS */
    label pcefl = 'PCE flag (<5% vs >= 5%)';
    pcefl = ifc(pce < 5, 'LT5', 'GTE5');

    /* DEMOGRAPHICS */
    length agegp $5;
    label agegp = "Under 50 years of age (40-49) or not (50+)";
    agegp = ifc(age < 50, '40-49', '50+');


    length gbdrace $30;
    label gbdrace = 'GBD + Race (if hi-income) or GBD only (if not)';
    gbdrace = '';
    if gbdgp1 ne '' then do;
        if gbdgp1 = 'High Income' then do;
            if racegp='White' then gbdrace=compress('High Income White');
            else if racegp='Black or African American' then gbdrace=compress('High Income Black');
            else gbdrace=compress('High Income Other');
        end;
        else gbdrace = compress(gbdgp1);
    end;

    /* CARDIOMETABOLIC */
    label cmdiabfl = 'Updated antidiabetes medication flag (Y/N)';
    cmdiabfl = ifc(_cmdiabfl='Y', 'Y', 'N');
    drop _cmdiabfl;

    label cmhypfl = 'Updated antihypertensive medication flag (Y/N)';
    cmhypfl = ifc(_cmhypfl='Y', 'Y', 'N');
    drop _cmhypfl;


    label diabfl = 'Updated pre-existing diabetes flag (Y/N)';
    diabfl = ifc(_diabfl='Y', 'Y', 'N');
    drop _diabfl;

    label egfrfl = 'EGFR flag (<90 vs >= 90)';
    length egfrfl $5;
    egfrfl= '';
    if egfr_eq21 ne . then do;
        egfrfl = ifc(egfr_eq21 < 90, 'LT90', 'GTE90');
    end; 

    label fhcvd2 = 'Updated family history of CVD flag (Y/N) - Unknown=missing';
    fhcvd2 = ifc(fhcvd='Unknown', '', fhcvd);

    label hdlfl = 'Screening HDL flag (Low (F<40; M<50) vs High)';
    length hdlfl $5;
    hdlfl = '';
    if shdl = . then hdlfl = '';
    else if shdl ne . and natalsex='Female' then do;
        if shdl < 40 then hdlfl='Low';
        else hdlfl='High';
    end;
    else if shdl ne . and natalsex='Male' then do;
        if shdl < 50 then hdlfl='Low';
        else hdlfl='High';
    end;

    label fglucfl = 'Fasting glucose flag (<100 vs >= 100)';
    if fgluc ne . then do;
        fglucfl = ifc(fgluc < 100, 'LT100', 'GTE100');
    end;
    else fglucfl = '';

    label hypfl = 'Updated pre-existing hypertension flag (Y/N)';
    hypfl = ifc(_hypfl='Y', 'Y', 'N');
    drop _hypfl;

    label ldlfl = 'Screening LDL flag (Low (<100) vs High)';
    ldlfl = ifc(sldlc < 100, 'LT100', 'GTE100');

    label metS2 = 'Metabolic syndrome (Not evaluable = missing)';
    metS2 = ifc(metS='Not evaluable', '', metS);

    label prediab2 = 'Updated prediabetes category (compressed)';
    length prediab2 $5;
    if      prediab='Diabetic (MH)'     then prediab2='D';
    else if prediab='Pre-diabetic'      then prediab2='PD';
    else if prediab='Not pre-diabetic'  then prediab2='NPD';
    else if prediab='Undetermined'      then prediab2='U';

    label trigfl = 'Screening triglycerides flag (<150 vs >= 150)';
    trigfl = ifc(strig < 150, 'LT150', 'GTE150');

    /* HIV1 */
    label abcfl = 'Updated abacavir flag (Y=current or former/N)';
    abcfl = ifc(_abcfl='Y', 'Y', 'N');
    drop _abcfl;

    label instifl = 'Current INSTI flag (Y=current/N)';
    instifl = ifc(find(acat2, 'INSTI'), 'Y', 'N');

    label nadirgp2 = 'Updated nadir CD4 grouping (Unknown=missing)';
    nadirgp2_temp = ifc(_nadirgp2='Unknown', '', _nadirgp2);
    if nadirgp2_temp ne '' then do;
        nadirgp2 = ifc(nadirgp2_temp='<200', 'LT200', 'GTE200');
    end;
    else nadirgp2='';
    drop _nadirgp2 nadirgp2_temp;

    label arvdurgp = 'Updated lifetime ART duration (Unknown=missing)';
    arvdurgp = ifc(_arvdurgp='Unknown', '', _arvdurgp);
    drop _arvdurgp;

    /* LIFESTYLE */
    label cigstgp = 'Updated smoking status (compressed)';
    cigstgp0 = tranwrd(_cigstgp, 'Smoked', '');
    cigstgp  = tranwrd(cigstgp0, 'Smoker', '');
    drop _cigstgp cigstgp0;
run;

proc print data=adsl_derived(obs=10);
    title 'Preview ADSL + derived variables';
run;

/*******************************************************************************
                                M061_ADTTMACE                                
*******************************************************************************/
/* 

    PATID       Participant identifier
    PARAM       Parameter (event type)
    CNSR        Censor (1=censored, 2=competing event)
    EVNT        Event  (1=event)
    ADT         Analysis date (event date)
    STARTDT     TTE origin date (entry date)

*/

title 'm061_adttmace: Event types (PARAM)';
proc freq data=derived.m061_adttmace;
    table param;
run;

data adttmace;
    set derived.m061_adttmace(where=(PARAM='First Primary MACE'));
    * CODE COMPETING EVENT AS CENSORING;
    if CNSR=2 then CNSR=1;
    keep PATID PARAMCD PARAM CNSR EVNT ADT AVAL STARTDT;
run;

title 'Preview adttmace';
proc print data=adttmace(obs=10);
run;

/*******************************************************************************
                                M081_ADTTAE_FINAL                                
*******************************************************************************/
/* 

    PATID       Participant identifier
    PARAM       Parameter (event type)
    CNSR        Censor (1=censored, 2=competing event)
    EVNT        Event  (1=event)
    ADT         Analysis date (event date)
    STARTDT     TTE origin date (entry date)

*/

title 'm081_adttae_final: Event types (PARAM)';
proc freq data=derived.m081_adttae_final;
    table paramcd*param/ list;
run;

data adttae;
    set derived.m081_adttae_final(where=(PARAM='1st Occurrence of Confirmed Diabetes'));
    * CODE COMPETING EVENT AS CENSORING;
    if CNSR=2 then CNSR=1;
    keep PATID PARAMCD PARAM CNSR EVNT ADT AVAL STARTDT;
run;

title 'Preview adttae';
proc print data=adttae(obs=10);
run;

/*******************************************************************************
                                M19_ADBIOMRKRS                               
*******************************************************************************/

proc freq data=derived.m19_adbiomrkrs;
    table avisit*paramcd*param/list;
    where avisit='Entry';
run;

* Extract hs-CRP at baseline only;
data adbiomrkrs;
    set derived.m19_adbiomrkrs(where=(PARAMCD='HSCRP001' and AVISIT='Entry' and ANL03FL='Y') rename=(AVAL=hscrp));
    keep PATID hsCRP;
    label hscrp='hs-CRP (mg/L)';
run;

title 'Preview adbiomrkrs'; 
proc print data=adbiomrkrs(obs=10);
run;

/*******************************************************************************
                                M02_RISKSCORE                              
*******************************************************************************/

* Extract 10-year PREVENT ASCVD score only;
data riskscore;
    set derived.m02_riskscore;
    keep PATID pvt_ascvd10 pvt_ascvd10g1 pvt_lowfl bool3;
    if pvt_ascvd10 ne . then do;
        label pvt_lowfl="PREVENT 10-year risk score < 3%";
        pvt_lowfl=ifc(pvt_ascvd10 <3, 'Low risk', 'Not low risk');
        label bool3='Low risk boolean';
        bool3=ifc(pvt_ascvd10 <3, 'TRUE', 'FALSE');
    end;
run;

title 'Preview riskscore'; 
proc print data=riskscore(obs=10);
run;


/*******************************************************************************
                                ADGRF                                
*******************************************************************************/

data derived.adgrf;
    merge adsl_derived
          adttmace(keep=patid cnsr evnt aval rename=(cnsr=cnsrm evnt=evntm aval=avalm)) 
          adttae  (keep=patid cnsr evnt aval rename=(cnsr=cnsrd evnt=evntd aval=avald))
          adbiomrkrs
          riskscore;
    by patid;
    label cnsrm = 'Censor indicator: Primary MACE (1=cnsr)'
          cnsrd = 'Censor indicator: Confirmed Diabetes (1=cnsr)'
          avalm = 'TTE: Primary MACE (years)'
          avald = 'TTE: Confirmed Diabetes (years)'
          evntm = 'Event indicator: Primary MACE (1=evnt)'
          evntd = 'Event indicator: Confirmed Diabetes (1=evnt)';
run;

%PGMSRC(INDATA = derived.adgrf,
        CLASS  = ADSL,
        LABEL  = Generalized Random Forest (GRF) Inputs,
        MISLAB = Y,
        VARORD = N,
        FILE   = Y,
        NUMOBS = 60);

%PERMSET(INDATA = derived.adgrf);

proc contents data=derived.adgrf;
run;

proc print data=derived.adgrf(obs=10);
run;