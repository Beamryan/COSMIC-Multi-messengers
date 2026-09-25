// ******************************************************************
//
//  MODULE  : gxProcess_COSMIC_2024.c
//  MODIFIED: 27 Sep 2024
//
//  AUTHORS : Shane L. Larson (Northwestern University)
//           s.larson@northwestern.edu
//
//  COMPILE COMMAND: gcc -lm gxProcess_COSMIC_2021.c
//
//  OTHER REQUIRED FILES:
//     processInputs.txt -- input file with your choices for processing
//     Sensitivity Curve -- 2 column CSV file of sensitivity curve
//     hilsBenderBackground.csv
//     hlrkBackground.csv
//  Right now the length of these files is hard-coded so arrays can be
//  dynamically allocated using nCURVE and nHBCURVE; improvement will
//  be to have the length read out of the file.
//
//  VERSIONS (incomplete record) ------------------------------------
//  2023 Jun 16: Added galactic component to column records
//
//  2022 May 23: Inserts paths for directories to output files.
//
//  2022 Jun 20: Added 2 processing outputs -- my files, and just the
//               standard COSMIC file output but separated by processing
//
//  2024 Sep 27: Changed input file to use cosmicRuntimeValues.txt which
//               is a collected file with all parameters for a given
//               COSMIC run, rather than two files.
//
//  DESCRIPTION
//  ----------------------------------------------------------------
//   This program looks through a processed COSMIC galaxy, and computes
//   the SNR of each source against a standard sensitivity curve, and
//   sorts the stars into {MONOCHROME, CHIRPING, CONFUSED} and writes
//   separate files. It also makes a RESOLVED file of MONO + CHIRPING.
//
//   INPUT GALXY FILE       : COSMIC 3.3, https://cosmic-popsynth.github.io
//   CSV PRE-PROCESSOR CODE : convert_MW_compact_binaries_to_csv_v2.py (Julie Malewicz, July 2020)
//   INPUT CSV FILE         : single common csv file for entire galaxy
//
//   OUTPUT CSV FILE(S)     : Binary parameter files sorted by LISA view of source:
//                            gxBasename_confused.csv
//                            gxBasename_monochromatic.csv
//                            gxBasename_chirping.csv
//                            gxBasename_resolved.csv
//                            gxBasename_runtimeData.txt
//
// ******************************************************************


//  ============ INCLUDE LIBRARIES ============

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>



//  ============ CONSTANTS DEFINED ============

#define MSun            1.989E30                   //  solar mass in kg
#define AU              149597900000.0             //  m in an AU
#define lyr             9.46052840488E15           //  m in a lyr
#define pc              3.08567818585E16           //  m in a pc
#define yr              31558149.7635456           // seconds per sidereal year

#define PI              3.14159265358979323846264338327950288  //  Why ask pi?
#define deg2rad         1.7453292519943295769236907684886e-2   //  radians per degree
#define rad2deg         57.295779513082320876798154814105      //  degrees per radian

#define G               6.67408E-11                //  Newton G, mks PDG
#define C               299792458.0                //  Speed of light, mks

#define MAX_LENGTH      2048                        // max string length


//  ============ FUNCTION PROTOTYPES ============

double GetSensitivity(int wdFLAG, double fgw, double fL[], double hfL[], double fHB[], double hfHB[],
                      int nCURVE, int nHBCURVE, int Tyr);


void parseTokenString(char *myString, char *myRead, char *delimeter);
double parseTokenDouble(char *myRead, char *delimeter);
int parseTokenInt(char *myRead, char *delimeter);
long parseTokenLong(char *myRead, char *delimeter);


void ErrorExit(char routine[], char errorMsg[]);

//  =========================================================
//    MAIN  ROUTINE
//  =========================================================

int main(void)
{
    //  ============= VARIABLE DECLARATIONS =============
    
    int ii, jj, kk, qq, rr;   //  looping indexes
    int chirpFLAG, wdFLAG;
    int nCOL, nCURVE, nHBCURVE;
    int lineNum;
    
    unsigned long int binNum, nCOUNT, nREAD, nBRIGHT, nMONOCHROME, nCHIRP;
    int nYRS;
    double fgw, SNR, bigSNR;
    double *fL, *hfL, *fHB, *hfHB;
    
    double tphys, mass1, mass2, kstar1, kstar2, evolType;
    double sepFinal, eccFinal, Porb, fgwPeak;
    double xGx, yGx, zGx, dist;
    double rad1, rad2, lum1, lum2, teff1, teff2, spin1, spin2, tbirth, tevGW;
    double magF1, magF2;
    
    double m1, m2, mc, dSI;
    double ho, SNRthresh, Tmp1, Tmp2, Tmp3;
    double m, b, hf, hfSRC, hfLISA;
    double Tobs, ALPHA, kappa;
    double deltaf, dFreq, fdot, fNow, fEnd, tFreq, hoFreq, hfLISAc, sumSNR;
    double xSun, ySun, zSun, xSrc, ySrc, zSrc, r;
    
    // --- date & time processing ---
    time_t startTime, endTime;
    struct tm *timeQuery;
    double elapsedTime;
    char timeBuffer1[128], timeBuffer2[128];
    
    //  --- File Handling ---
    
    FILE *myFile, *inLISA, *gxParams, *gxMonochrome, *gxChirp, *gxConfused, *gxResolved;
    FILE *sllMonochrome, *sllConfused, *sllResolved, *sllChirp, *sllNS, *gxNS;
    FILE *sllNSbright, *gxNSbright;
    FILE *gxRunData;
    char basename[255], runFilename[255], wdFilename[255], sensitivityFilename[255];
    char galaxyFilename[255], confusedFilename[255], monochromeFilename[255], chirpFilename[255], resolvedFilename[255], nsFilename[255],nsBrightFilename[255];
    char confusedSLL[255], monochromeSLL[255], chirpSLL[255], resolvedSLL[255],
        nsSLL[255],nsBrightSLL[255];

    char csvOutputLine[1024], cosmicOutputLine[1024];

    char errorMsg[255];
    
    char tmpRead[MAX_LENGTH + 1];
    char *tokenParse, *myLine;  // strings for parsing tokens (elements) from a CSV line
    char tString;
    char gxComp[128];
    
    
    //  =============== START ROUTINE HERE ===============
    
    printf("Welcome to GALACTIC LISA PROCESSOR DELUXE.\n\n");
    
    //  -------- initialization parameters -------------------------------------------
    
    nREAD = 0;
    nCOUNT = 0;
    nBRIGHT = 0;
    nMONOCHROME = 0;
    nCHIRP = 0;
    
    bigSNR = 0.0;
    
    
    //  -----------------------------------------------------------------
    //  --- INPUT PARAMETERS READ FROM FILE BEFORE A RUN ----------------
    //  -----------------------------------------------------------------
    
    // FUTURE IMPROVEMENT: take processing filename as command line input
    if ((myFile = fopen("../buildFiles/cosmicRuntimeValues.txt","r")) == 0)
        ErrorExit("main","File   cosmicRuntimeValues.txt  NOT FOUND");
    
    lineNum = 0;
    while (fgets(tmpRead,MAX_LENGTH + 1,myFile))  //  as long as not at EOF
    {
        lineNum++;
        myLine = tmpRead;                        // set up string for token parsing
        
        if (lineNum == 6) parseTokenString(basename,tmpRead,"=");
        if (lineNum == 51) parseTokenString(sensitivityFilename,tmpRead,"=");
        
        if (lineNum == 52) nYRS = parseTokenLong(tmpRead,"=");         // Tobs
        if (lineNum == 53) SNRthresh = parseTokenDouble(tmpRead,"=");  // THRESHOLD SNR for detection
        if (lineNum == 54) ALPHA = parseTokenDouble(tmpRead,"=");      // CHIRP THRESHOLD -- number of bins for Tobs which is detectable
        if (lineNum == 55) wdFLAG = parseTokenInt(tmpRead,"=");        // 0 = NO WD;   1 = Cornish Robson;    2 = HLRK;    3 = Hils Bender

    }
    fclose(myFile);

    Tobs = (double)nYRS*yr;    //  LISA Observing Time in SECONDS -- sets bin width
    
    
    // -------------------------------------------------------------------------
    //  --- Get LISA Sensitivity Curve -----------------------------------------
    // -------------------------------------------------------------------------
    
    if ((inLISA = fopen(sensitivityFilename,"r")) == 0)       // sensitivity curve file
    {
        sprintf(errorMsg,"Sensitivity File   %s   NOT FOUND",sensitivityFilename);
        ErrorExit("main",errorMsg);
    }
    
    nCOL = 2;                                // number of columns in sensitivity & wdBkgnd files
    ii = 0;                                  // data point counter for LISA arrays
    tString = '#';                           // header line intro character for search
    
    // --- count lines in sensitivity curve file, skipping over header lines
    nCURVE = 0;
    while (fgets(tmpRead,MAX_LENGTH + 1,inLISA))  //  as long as not at EOF
    {
        myLine = tmpRead;                        // set up string for token parsing
        if (strchr(tmpRead,tString) == NULL)      // check for header line character
        {
            nCURVE++;                           // increment data point counter
        }
    }
    rewind(inLISA);
    
    // -- allocate arrays for sensitivity ----
    fL = (double *) malloc(nCURVE*sizeof(double));
    hfL = (double *) malloc(nCURVE*sizeof(double));
    
    // --- read in sensitivity curve file, skipping over header lines
    while (fgets(tmpRead,MAX_LENGTH + 1,inLISA))  //  as long as not at EOF
    {
        myLine = tmpRead;                         // set up string for token parsing
        
        if (strchr(tmpRead,tString) == NULL)      // check for header line character
        {
            for (qq = 0; qq < nCOL; qq++)         // loop over number of columns in the csv file
            {
                tokenParse = strsep(&myLine,","); // look for next token -- separator is ","
                
                if (qq == 0) fL[ii] = atof(tokenParse);     // first column is frequency
                if (qq == 1) hfL[ii] = atof(tokenParse);    // second column is spectral amplitude
            }
            
            // printf("%d\n",ii);                    // echo
            // fflush(stdout);                       // echo
            ii++;                              // increment data point counter
        }
    }
    
    fclose(inLISA);
    
    printf("Got sensitivity...\n");
    fflush(stdout);
    
    //  -----------------------------------------------------------------
    //  --- LISA CONFUSION FILES (if needed) ----------------------------
    //  -----------------------------------------------------------------
    
    //  --- Hils Bender Galactic Background ---
    if ((wdFLAG == 2) || (wdFLAG == 3))
    {
        if (wdFLAG == 2)
            strcpy(wdFilename,"../lisaFiles/hlrkBackground.csv");         // HLRK fit to HilsBender
        else
            strcpy(wdFilename,"../lisaFiles/hilsBenderBackground.csv");    // Hils Bender
        
        
        if ((inLISA = fopen(wdFilename,"r")) == 0)       // file to read in pre-defined background
        {
            sprintf(errorMsg,"White Dwarf Noise File   %s   NOT FOUND",wdFilename);
            ErrorExit("main",errorMsg);
        }
        

        
        // --- count lines in confusion curve file, skipping over header lines
        nHBCURVE = 0;
        while (fgets(tmpRead,MAX_LENGTH + 1,inLISA))  //  as long as not at EOF
        {
            myLine = tmpRead;                        // set up string for token parsing
            if (strchr(tmpRead,tString) == NULL)      // check for header line character
            {
                nHBCURVE++;                           // increment data point counter
            }
        }
        rewind(inLISA);

        // -- allocate arrays for confusion curve ----
        fHB = (double *) malloc(nHBCURVE*sizeof(double));
        hfHB = (double *) malloc(nHBCURVE*sizeof(double));
        
        
        // ---------------------------------------------------------------
        // --- BACKGROUND CURVE FILE read, skipping over header lines
        // ---------------------------------------------------------------
        while (fgets(tmpRead,MAX_LENGTH + 1,inLISA))  //  as long as not at EOF
        {
            myLine = tmpRead;                         // set up string for token parsing
            
            if (strchr(tmpRead,tString) == NULL)      // check for header line character
            {
                for (qq = 0; qq < nCOL; qq++)         // loop over number of columns in the csv file
                {
                    tokenParse = strsep(&myLine,","); // look for next token -- separator is ","
                    
                    if (qq == 0) fHB[ii] = atof(tokenParse);     // first column is frequency in Hz
                    if (qq == 1) hfHB[ii] = atof(tokenParse);    // second column is spectral amplitude in 1/rtHz
                }
                
                // printf("%d\n",ii);                    // echo
                // fflush(stdout);                       // echo
                ii++;                              // increment data point counter
            }
        }
        
        fclose(inLISA);
    } // END --->  HLRK or HilsBender is chosen
    
    
    
    // -------------------------------------------------------------------------
    // ---- output files setup -------------------------------------
    // -------------------------------------------------------------------------

    // --- output, COSMIC CSV record -------------------------------
    sprintf(runFilename,"../errLogs/%s_runtimeData.log",basename);
    sprintf(confusedFilename,"../csv/%s_confused.csv",basename);
    sprintf(monochromeFilename,"../csv/%s_monochrome.csv",basename);
    sprintf(chirpFilename,"../csv/%s_chirping.csv",basename);
    sprintf(resolvedFilename,"../csv/%s_resolved.csv",basename);
    sprintf(galaxyFilename,"../csv/%s_full.csv",basename);
    sprintf(nsFilename,"../csv/%s_nsAll.csv",basename);
    sprintf(nsBrightFilename,"../csv/%s_nsBright.csv",basename);

    // --- output, SLL CSV record ----------------------------------
    sprintf(confusedSLL,"../csv/sll_%s_confused.csv",basename);
    sprintf(monochromeSLL,"../csv/sll_%s_monochrome.csv",basename);
    sprintf(chirpSLL,"../csv/sll_%s_chirping.csv",basename);
    sprintf(resolvedSLL,"../csv/sll_%s_resolved.csv",basename);
    sprintf(nsSLL,"../csv/sll_%s_nsAll.csv",basename);
    sprintf(nsBrightSLL,"../csv/sll_%s_nsBright.csv",basename);

    
    
    // -- open main galaxy up -----------------
    if ((gxParams = fopen(galaxyFilename,"r")) == 0)       // file to read in galaxy
    {
        sprintf(errorMsg,"GALAXY Data File   %s   NOT FOUND",galaxyFilename);
        ErrorExit("main",errorMsg);
    }
    
    // --- output files open up ----
    gxRunData = fopen(runFilename,"w");
    gxConfused = fopen(confusedFilename,"w");
    gxMonochrome = fopen(monochromeFilename,"w");
    gxResolved = fopen(resolvedFilename,"w");
    gxChirp = fopen(chirpFilename,"w");
    gxNS = fopen(nsFilename,"w");
    gxNSbright = fopen(nsBrightFilename,"w");

    sllConfused = fopen(confusedSLL,"w");
    sllMonochrome = fopen(monochromeSLL,"w");
    sllResolved = fopen(resolvedSLL,"w");
    sllChirp = fopen(chirpSLL,"w");
    sllNS = fopen(nsSLL,"w");
    sllNSbright = fopen(nsBrightSLL,"w");

    
    // --- output meta data to runtime file ------------
    // -------------------------------------------------
    
    // --- get the time of the run ---
    time(&startTime);                                     // get time from the system
    timeQuery = localtime(&startTime);                    // process time
    strftime(timeBuffer1,128,"%x - %I:%M%p",timeQuery);  // human readable string
                                                         // printf("Time is: %s\n",timeBuffer);
   
    printf("\n ========================================== \n");
    printf("\n Start Time   : %s\n",timeBuffer1);

    // --- output headers ------
    fprintf(gxRunData,"# CODE: gxProcess_COSMIC_2022.c    sll, June 2023\n");
    fprintf(gxRunData,"#\n");
    fprintf(gxRunData,"# Galaxy Input File: %s.csv\n",basename);
    fprintf(gxRunData,"#\n");
    fprintf(gxRunData,"# Output File Columns are COMMA DELIMITED. Headers and Units are: \n");
    fprintf(gxRunData,"# binNum, gxComponent, mass1(mSun), mass2(mSun),lum1(Lsun),lum2(Lsun),Teff1(K),Teff2(K),spin1(rad/yr),spin2(rad/yr),Porb(s),xGx(pc),yGx(pc),zGx(pc),D(pc),r(pc),fgw(Hz),fdot(Hz/s),ho,hf(per rtHz),SNR\n");

    
    // --- print headers to files ---
    fprintf(gxConfused,"#bin_num,tphys,mass_1(Msun),mass_2(Msun),kstar_1,kstar_2,evol_type,rad_1(Rsun),rad_2(Rsun),lum_1(Lsun),lum_2(Lsun),teff_1(K),teff_2(K),omega_spin_1(yr^-1),omega_spin_2(yr^-1),B_1(Gauss),B_2(Gauss),t_birth,t_evol_GW,sep_final(Rsun),ecc_final,porb_final(day),xGx(kpc),yGx(kpc),zGx(kpc),dist(kpc),gxComponent\n");
    fprintf(gxMonochrome,"#bin_num,tphys,mass_1(Msun),mass_2(Msun),kstar_1,kstar_2,evol_type,rad_1(Rsun),rad_2(Rsun),lum_1(Lsun),lum_2(Lsun),teff_1(K),teff_2(K),omega_spin_1(yr^-1),omega_spin_2(yr^-1),B_1(Gauss),B_2(Gauss),t_birth,t_evol_GW,sep_final(Rsun),ecc_final,porb_final(day),xGx(kpc),yGx(kpc),zGx(kpc),dist(kpc),gxComponent\n");
    fprintf(gxChirp,"#bin_num,tphys,mass_1(Msun),mass_2(Msun),kstar_1,kstar_2,evol_type,rad_1(Rsun),rad_2(Rsun),lum_1(Lsun),lum_2(Lsun),teff_1(K),teff_2(K),omega_spin_1(yr^-1),omega_spin_2(yr^-1),B_1(Gauss),B_2(Gauss),t_birth,t_evol_GW,sep_final(Rsun),ecc_final,porb_final(day),xGx(kpc),yGx(kpc),zGx(kpc),dist(kpc),gxComponent\n");
    fprintf(gxResolved,"#bin_num,tphys,mass_1(Msun),mass_2(Msun),kstar_1,kstar_2,evol_type,rad_1(Rsun),rad_2(Rsun),lum_1(Lsun),lum_2(Lsun),teff_1(K),teff_2(K),omega_spin_1(yr^-1),omega_spin_2(yr^-1),B_1(Gauss),B_2(Gauss),t_birth,t_evol_GW,sep_final(Rsun),ecc_final,porb_final(day),xGx(kpc),yGx(kpc),zGx(kpc),dist(kpc),gxComponent\n");
    fprintf(gxNS,"#bin_num,tphys,mass_1(Msun),mass_2(Msun),kstar_1,kstar_2,evol_type,rad_1(Rsun),rad_2(Rsun),lum_1(Lsun),lum_2(Lsun),teff_1(K),teff_2(K),omega_spin_1(yr^-1),omega_spin_2(yr^-1),B_1(Gauss),B_2(Gauss),t_birth,t_evol_GW,sep_final(Rsun),ecc_final,porb_final(day),xGx(kpc),yGx(kpc),zGx(kpc),dist(kpc),gxComponent\n");
    fprintf(gxNSbright,"#bin_num,tphys,mass_1(Msun),mass_2(Msun),kstar_1,kstar_2,evol_type,rad_1(Rsun),rad_2(Rsun),lum_1(Lsun),lum_2(Lsun),teff_1(K),teff_2(K),omega_spin_1(yr^-1),omega_spin_2(yr^-1),B_1(Gauss),B_2(Gauss),t_birth,t_evol_GW,sep_final(Rsun),ecc_final,porb_final(day),xGx(kpc),yGx(kpc),zGx(kpc),dist(kpc),gxComponent\n");


    fprintf(sllConfused,"# binNum,gxComponent,kstar1,kstar2,mass1(mSun),mass2(mSun),lum1(Lsun),lum2(Lsun),Teff1(K),Teff2(K),rad_1(Rsun),rad_2(Rsun),B1(Gauss),B2(Gauss),spin1(rad/yr),spin2(rad/yr),Porb(s),xGx(pc),yGx(pc),zGx(pc),D(pc),fgw(Hz),fdot(Hz/s),ho,hf(per rtHz),SNR\n");
    fprintf(sllMonochrome,"# binNum,gxComponent,kstar1,kstar2,mass1(mSun),mass2(mSun),lum1(Lsun),lum2(Lsun),Teff1(K),Teff2(K),rad_1(Rsun),rad_2(Rsun),B1(Gauss),B2(Gauss),spin1(rad/yr),spin2(rad/yr),Porb(s),xGx(pc),yGx(pc),zGx(pc),D(pc),fgw(Hz),fdot(Hz/s),ho,hf(per rtHz),SNR\n");
    fprintf(sllChirp,"# binNum,gxComponent,kstar1,kstar2,mass1(mSun),mass2(mSun),lum1(Lsun),lum2(Lsun),Teff1(K),Teff2(K),rad_1(Rsun),rad_2(Rsun),B1(Gauss),B2(Gauss),spin1(rad/yr),spin2(rad/yr),Porb(s),xGx(pc),yGx(pc),zGx(pc),D(pc),fgw(Hz),fdot(Hz/s),ho,hf(per rtHz),SNR\n");
    fprintf(sllResolved,"# binNum,gxComponent,kstar1,kstar2,mass1(mSun),mass2(mSun),lum1(Lsun),lum2(Lsun),Teff1(K),Teff2(K),rad_1(Rsun),rad_2(Rsun),B1(Gauss),B2(Gauss),spin1(rad/yr),spin2(rad/yr),Porb(s),xGx(pc),yGx(pc),zGx(pc),D(pc),fgw(Hz),fdot(Hz/s),ho,hf(per rtHz),SNR\n");
    fprintf(sllNS,"# binNum,gxComponent,kstar1,kstar2,mass1(mSun),mass2(mSun),lum1(Lsun),lum2(Lsun),Teff1(K),Teff2(K),rad_1(Rsun),rad_2(Rsun),B1(Gauss),B2(Gauss),spin1(rad/yr),spin2(rad/yr),Porb(s),xGx(pc),yGx(pc),zGx(pc),D(pc),fgw(Hz),fdot(Hz/s),ho,hf(per rtHz),SNR\n");
    fprintf(sllNSbright,"# binNum,gxComponent,kstar1,kstar2,mass1(mSun),mass2(mSun),lum1(Lsun),lum2(Lsun),Teff1(K),Teff2(K),rad_1(Rsun),rad_2(Rsun),B1(Gauss),B2(Gauss),spin1(rad/yr),spin2(rad/yr),Porb(s),xGx(pc),yGx(pc),zGx(pc),D(pc),fgw(Hz),fdot(Hz/s),ho,hf(per rtHz),SNR\n");


    // -------------------------------------------------------------------------
    // ---- PROCESS GALAXY FILE  -------------------------------------
    // -------------------------------------------------------------------------

    // This is the token parsing for the original COSMIC 3.3 Malewicz file
    // we were using in 2020.
    //
    // if (qq == 0) binNum = atol(tokenParse);     // col  0 = binnum
    // if (qq == 1) tphys = atof(tokenParse);      // col  1 = tphys (Myr)
    // if (qq == 2) mass1 = atof(tokenParse);      // col  2 = mass 1 (Msun)
    // if (qq == 3) mass2 = atof(tokenParse);      // col  3 = mass 2 (Msun)
    // if (qq == 4) kstar1 = atof(tokenParse);     // col  4 = kstar
    // if (qq == 5) kstar2 = atof(tokenParse);     // col  5 = kstar
    // if (qq == 6) sepFinal = atof(tokenParse);   // col  6 = separation (Rsun)
    // if (qq == 7) eccFinal = atof(tokenParse);   // col  7 = eccentricity
    // if (qq == 8) Porb = atof(tokenParse);       // col  8 = orbital period (days)
    // if (qq == 9) fgwPeak = atof(tokenParse);    // col  9 = Peak GW Freq (Hz)
    // if (qq == 10) xGx = atof(tokenParse);       // col  10 = galactocentric (kpc)
    // if (qq == 11) yGx = atof(tokenParse);       // col  11 = galactocentric (kpc)
    // if (qq == 12) zGx = atof(tokenParse);       // col  12 = galactocentric (kpc)
    // if (qq == 13) dist = atof(tokenParse);      // col  13 = Distance from Sun (kpc)


    // --- galaxy file read initilization data
    nCOL = 28;                                // number of columns in cosmic galaxy CSV file
    tString = '#';                            // header line intro character for search
    
    // --- read in galaxy data file, echo header lines to process files
    while (fgets(tmpRead,MAX_LENGTH + 1,gxParams))  //  as long as not at EOF
    {
        myLine = tmpRead;                         // set up string for token parsing
        strcpy(cosmicOutputLine,tmpRead);         // copy the COSMIC csv line for output
        
        //  FUTURE IMPROVEMENT: : echo any header lines to the processed output
        //  files, so the galaxy data propagates with the processed files
        
        if (strchr(tmpRead,tString) == NULL)      // check for header line character
        {
            for (qq = 0; qq < nCOL; qq++)         // loop over number of columns in the csv file
            {
                tokenParse = strsep(&myLine,","); // look for next token -- separator is ","
                
                if (qq == 0) binNum = atol(tokenParse);     // col  0 = binnum
                if (qq == 1) tphys = atof(tokenParse);      // col  1 = tphys (Myr)
                if (qq == 2) mass1 = atof(tokenParse);      // col  2 = mass 1 (Msun)
                if (qq == 3) mass2 = atof(tokenParse);      // col  3 = mass 2 (Msun)
                if (qq == 4) kstar1 = atof(tokenParse);     // col  4 = kstar
                if (qq == 5) kstar2 = atof(tokenParse);     // col  5 = kstar
                if (qq == 6) evolType = atof(tokenParse);   // col  6 = evol_type

                if (qq == 7) rad1 = atof(tokenParse);       // col  7 = radius 1 (Rsun)
                if (qq == 8) rad2 = atof(tokenParse);       // col  8 = radius 2 (Rsun)

                if (qq == 9) lum1 = atof(tokenParse);       // col  7 = luminosity 1 (Lsun)
                if (qq == 10) lum2 = atof(tokenParse);       // col  8 = luminosity 2 (Lsun)
                if (qq == 11) teff1 = atof(tokenParse);      // col  9 = Teff 1 (K)
                if (qq == 12) teff2 = atof(tokenParse);      // col  10 = Teff 2 (K)
                if (qq == 13) spin1 = atof(tokenParse);     // col  11 = omega spin 1 (rad/yr)
                if (qq == 14) spin2 = atof(tokenParse);     // col  12 = omega spin 2 (rad/yr)
                if (qq == 15) magF1 = atof(tokenParse);     // col  13 = B field 1 (Gauss)
                if (qq == 16) magF2 = atof(tokenParse);     // col  14 = B field 1 (Gauss)
                if (qq == 17) tbirth = atof(tokenParse);    // col  15 = tbirth (Myr)
                if (qq == 18) tevGW = atof(tokenParse);     // col  16 = tevolveGW (Myr)

                if (qq == 19) sepFinal = atof(tokenParse);  // col  17 = separation (Rsun)
                if (qq == 20) eccFinal = atof(tokenParse);  // col  18 = eccentricity
                if (qq == 21) Porb = atof(tokenParse);      // col  19 = orbital period (days)
                if (qq == 22) xGx = atof(tokenParse);       // col  20 = galactocentric (kpc)
                if (qq == 23) yGx = atof(tokenParse);       // col  22 = galactocentric (kpc)
                if (qq == 24) zGx = atof(tokenParse);       // col  22 = galactocentric (kpc)
                if (qq == 25) dist = atof(tokenParse);      // col  23 = Distance from Sun (kpc)
                if (qq == 26) strcpy(gxComp,tokenParse);    // col  24 = gx Component: ThinDisk,ThickDisk,Bulge
                if (qq == 27) gxComp[strcspn(gxComp, "\r\n")] = 0;  //removes trailing carriage return

            }

            nREAD++;    // increase count of read binaries
            nCOUNT++;   // increment successul read counter --> for SCREEN PROGRESS OUTPUT at end of this loop
            
            // ---- do the GW calculation from the source parameters ------------------
            // ------------------------------------------------------------------------
            
            //  convert to heliocentric position ----------------------------
            //   SUN'S COORDINATES:
            //     xSun = 7860 pc  [ Boehle et al, ApJ 830, 17 (2016) ]
            //     ySun = 0 pc     [ assumption of coordinate system ]
            //     zSun = 27 pc    [ Chen et al, ApJ 553, 184 (2001) ]
            
            xSun = 7860.0;    // I think this is opposite astropy, which assumes Sun on -x axis
            ySun = 0.0;
            zSun = 27.0;
            
            // convert to coordinates centered on the Sun (LISA)
            //   the *1000 converts from kpc (kmb reports) to pc
            
            xSrc = xGx*1000.0 - xSun;
            ySrc = yGx*1000.0 - ySun;
            zSrc = zGx*1000.0 - zSun;
            
            r = sqrt(xSrc*xSrc + ySrc*ySrc + zSrc*zSrc)*pc;    // distance to soruce in meters
            dSI = dist*1000.0*pc;                               // cosmic distance in meters
            r = dSI;                                           // use COSMIC distance
            
            // parameter calculations ----
            m1 = mass1*(MSun*G/(C*C));                       //  convert to meters
            m2 = mass2*(MSun*G/(C*C));                       //  convert to meters
            mc = (pow(m1*m2,3.0/5.0))/(pow(m1+m2,1.0/5.0));   // chrip mass in meters
            
            fgw = 2.0/(Porb*86400.0);                         //  GW frequency in Hz, CIRCULAR ORBITS
            
            // pocket formulae -----------
            ho = 4.0*(mc/r)*pow(PI*fgw*mc/C,2.0/3.0);               //  scaling amplitude
            hfSRC = sqrt(Tobs)*ho;                                 //  spectral amplitude of the source
            fdot = (96.0/5.0)*fgw*pow(PI*fgw*mc/C,8.0/3.0)/(mc/C);   //  quadrupole chirp in Hz/s
            
            // compute frequency at end of Tobs
            kappa = (5.0/256.0)*pow(mc/C,-5.0/3.0)*pow(PI,-8.0/3.0);
            
            Tmp1 = pow(fgw,-8.0/3.0) - Tobs/kappa;
            fEnd = pow(Tmp1,-3.0/8.0);             // find this by integrating quadrupole chirp
            
            // check if chirping ---------
            deltaf = ALPHA/Tobs;     // required shift in frequency for chirp detection
            
            if ( (fEnd - fgw) >= deltaf)
            {
                chirpFLAG = 1;  //  Chirping Binary!
            }
            else
            {
                chirpFLAG = 0;  //  Monochrome Binary!
            }
            
            // --- Get SNR ------------------------------------------------------------
            // ------------------------------------------------------------------------
            
            if ((fgw < fL[0]) || (fgw > fL[nCURVE-1]))
            {
                SNR = 0.0;    // if source does not overlap sensitivity, set SNR = 0
            }
            else
            {
                if (chirpFLAG == 0)
                {
                    hfLISA = GetSensitivity(wdFLAG, fgw, fL, hfL, fHB, hfHB, nCURVE, nHBCURVE, nYRS);
                    
                    SNR = hfSRC/hfLISA;          // simple SNR in line limit
                }
                else if (chirpFLAG == 1)
                {
                    
                    // Sum up all the contributions at each freq to get the SNR
                    sumSNR = 0.0;
                    fNow = fgw;                   // starting (input) frequency
                    dFreq = (fEnd - fgw)/100.0;   // divide up the total chirp range into 100 small bits for quick SNR calc
                    
                    for (kk = 0; kk < 100; kk++)    // loop over frequencies covered by chirping in the 100 small bits
                    {
                        tFreq = kappa*(pow(fNow,-8.0/3.0) - pow(fNow + dFreq,-8.0/3.0)); // time between frequencies
                        
                        hoFreq = 4.0*(mc/r)*pow(PI*fNow*mc/C,2.0/3.0);               //  scaling amplitude at the current freqeuncy
                        hfLISAc = GetSensitivity(wdFLAG, fNow, fL, hfL, fHB, hfHB, nCURVE, nHBCURVE, nYRS); // LISA value at current freqeuncy
                        
                        sumSNR += (hoFreq*hoFreq)*tFreq/(hfLISAc*hfLISAc);  // add up the contributions from each bin --> note tFreq makes numerator specAmp
                        
                        fNow += dFreq;  // update frequency to next little bit
                    }
                    
                    SNR = sqrt(sumSNR);    // chirping SNR value for the current binary
                }
            }

            
            // --- Process star based on SNR ------------------------------------------
            // ------------------------------------------------------------------------
            //  output parameters to file ------------------------
            //    columns in file are:
            //    col  1 = binnum
            //    col  2 = m1 (mSun)
            //    col  3 = m2 (mSun)
            //    col  4 = lum1 (Lsun)
            //    col  5 = lum2 (Lsun)
            //    col  6 = teff1 (K)
            //    col  7 = teff2 (K)
            //    col  8 = omegaSpin1 (rad/yr)
            //    col  9 = omegaSpin1 (rad/yr)
            //    col 10 = Porb (sec)
            //    col 11 = xGx (pc)
            //    col 12 = yGx (pc)
            //    col 13 = zGx (pc)
            //    col 14 = dist (pc)
            //    col 15 = fgw (Hz)
            //    col 16 = fdot (Hz/s)
            //    col 17 = ho
            //    col 18 = hf (per rtHz)
            //    col 19 = SNR
            // ---------------------------------------------------

            if (SNR > bigSNR) bigSNR = SNR; //  store biggest SNR
            
            // the output to all the files is THE SAME, just directed based on values, so make the line ONCE, to make output
            // and maintenance simpler -- use SPRINTF() -- there are 20 output columns
            
            sprintf(csvOutputLine,"%ld,%s,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g,%15.7g",binNum,gxComp,kstar1,kstar2,mass1,mass2,lum1,lum2,teff1,teff2,rad1,rad2,magF1,magF2,spin1,spin2,Porb*86400.0,xGx*1000.0,yGx*1000.0,zGx*1000.0,dist*1000.0,fgw,fdot,ho,hfSRC,SNR);
            
            //  If a source beats the request threshold, then store it and ouput data
            if (SNR > SNRthresh)
            {
                nBRIGHT++;  // update number of stars that are bright
                
                if (chirpFLAG == 0)
                {
                    fprintf(gxMonochrome,"%s",cosmicOutputLine);
                    fprintf(gxResolved,"%s",cosmicOutputLine);
                    fprintf(sllMonochrome,"%s\n",csvOutputLine);
                    fprintf(sllResolved,"%s\n",csvOutputLine);
                    nMONOCHROME++;  // increase count of monochrome binaries
                }
                
                if (chirpFLAG == 1)
                {
                    fprintf(gxChirp,"%s",cosmicOutputLine);
                    fprintf(gxResolved,"%s",cosmicOutputLine);
                    fprintf(sllChirp,"%s\n",csvOutputLine);
                    fprintf(sllResolved,"%s\n",csvOutputLine);
                    nCHIRP++;  // increase count of chirping binaries
                }
                
                // this dumps bright NS to file
                if (kstar1 == 13)
                {
                    fprintf(gxNSbright,"%s",cosmicOutputLine);
                    fprintf(sllNSbright,"%s\n",csvOutputLine);
                }
            }
            else
            {
                fprintf(gxConfused,"%s",cosmicOutputLine);
                fprintf(sllConfused,"%s\n",csvOutputLine);
            }
  
            // this dumps all NS to file
/*            if (kstar1 == 13)
            {
                fprintf(gxNS,"%s",cosmicOutputLine);
                fprintf(sllNS,"%s\n",csvOutputLine);
            }
*/
        }  // END --> processing non-header lines

        if (nCOUNT == 100000)
        {
            nCOUNT = 0;
            printf(".");
            fflush(stdout);
        }
        
        
    }   // END WHILE --> galaxy parameter file read
    
    fclose(gxParams);
    fclose(gxMonochrome);
    fclose(gxChirp);
    fclose(gxConfused);
    fclose(gxResolved);
    fclose(gxNS);
    fclose(gxNSbright);
    
    fclose(sllMonochrome);
    fclose(sllChirp);
    fclose(sllConfused);
    fclose(sllResolved);
    fclose(sllNS);
    fclose(sllNSbright);
    
    // --- final runtime data to file ----
    // --- get the end time of the run ---
    time(&endTime);                                     // get time from the system
    timeQuery = localtime(&endTime);                    // process time
    strftime(timeBuffer2,128,"%x - %I:%M%p",timeQuery);  // human readable string
    
    elapsedTime = difftime(endTime,startTime);
    
    fprintf(gxRunData,"\n Start Time   : %s\n",timeBuffer1);
    fprintf(gxRunData," End Time     : %s\n",timeBuffer2);
    fprintf(gxRunData," Process Time : %lf minutes\n",elapsedTime/60.0);
    
    
    fprintf(gxRunData,"\nSNR Threshold   : %lf\n",SNRthresh);
    fprintf(gxRunData,"\nTOTAL BINARIES    = %ld\n",nREAD);
    fprintf(gxRunData,"TOTAL RESOLVED    = %ld\n",nBRIGHT);
    fprintf(gxRunData,"TOTAL MONOCHROME  = %ld\n",nMONOCHROME);
    fprintf(gxRunData,"TOTAL CHIRPING    = %ld\n",nCHIRP);

    fclose(gxRunData);
    
    printf("\n\nTOTAL BINARIES    = %ld\n",nREAD);
    printf("ELAPSED TIME      = %lf minutes\n",elapsedTime/60.0);
    printf("TOTAL RESOLVED    = %ld\n",nBRIGHT);
    printf("TOTAL MONOCHROME  = %ld\n",nMONOCHROME);
    printf("TOTAL CHIRPING    = %ld\n",nCHIRP);
    printf("\n ========================================== \n");
    printf("\n\nALL DONE!  Whoo hoo!\n");
    
    return 0;
}


// ***************************************************************************
//
//  FUNCTION:  GetSensitivity
//
// ***************************************************************************

double GetSensitivity(int wdFLAG, double fgw, double fL[], double hfL[], double fHB[], double hfHB[],
                      int nCURVE, int nHBCURVE, int Tyr)
{
    // =============== VARIABLES ===============
    
    int jj;
    
    double m, b;
    double hfi, hfiL, hfiHB;     // interpolated spectral amplitude values
    
    double Sc, Ab, aCR, bCR, kCR, gCR, fk;  // Cornish Robson wd noise fit parameters
    
    double Tmp1, Tmp2;
    
    // =============== START ROUTINE HERE ===============
    
    // -----------------------------------------------------------------------
    //  --------- COMPUTE THE LISA CURVE VALUE FOR SNR ESTIMATE --------------
    
    //  look for frequency of interest in LISA curve
    if ((fgw > fL[0]) && (fgw < fL[nCURVE-1]))  // make sure the src is in the band covered by the sensitivity curve
    {
        jj = 0;
        while (fL[jj] <= fgw) jj++;
        
        //  straight line fit on log-log graph
        
        m = (log10(hfL[jj]/hfL[jj-1]))/(log10(fL[jj]/fL[jj-1])); //  Slope btwn LISA pts
        b = log10(hfL[jj]) - m*log10(fL[jj]);    //  Intercept for fit btwn LISA pts
        
        hfiL = pow(10.0,m*log10(fgw) + b);         //  Interpolated LISA value
        
        hfi = hfiL;    //  set value we use to LISA value, in case there is no WD contribution
    }
    
    
    // -----------------------------------------------------------------------
    //  --------- WHITE DWARF BACKGROUND VALUE FOR SNR ESTIMATE --------------
    //     WD FLAG -->  0 = NO WD;   1 = Cornish Robson;    2 = HLRK;    3 = Hils Bender
    
    if (wdFLAG == 1)
    {
        //  --- Cornish Robson Galactic Background ---
        Ab = 1.80e-44;
        
        if (Tyr == 1)
        {
            // --- 1 year parameters ---
            aCR = 0.171;
            bCR = 292.0;
            kCR = 1020.0;
            gCR = 1680.0;
            fk = 0.00215;
            
        }
        if (Tyr == 2)
        {
            // --- 2 year parameters ---
            aCR = 0.165;
            bCR = 299.0;
            kCR = 611.0;
            gCR = 1340.0;
            fk = 0.00173;
        }
        if (Tyr == 4)
        {
            // --- 4 year parameters ---
            aCR = 0.138;
            bCR = -221.0;
            kCR = 521.0;
            gCR = 1680.0;
            fk = 0.00113;
        }
        if ((Tyr != 1) && (Tyr != 2) && (Tyr !=4))
            ErrorExit("GetSensitivity","YRS do not match CornishRobson Choices");
        
        // --- Eq.3  in arxiv/1703.09858
        Tmp1 = pow(fgw,-7.0/3.0)*(1.0 + tanh(gCR*(fk - fgw)));
        Tmp2 = exp(-1.0*pow(fgw,aCR) + bCR*fgw*sin(kCR*fgw));
        Sc = Ab*Tmp1*Tmp2;
        
        hfiHB = sqrt(Sc);
        hfi = sqrt(hfiL*hfiL + hfiHB*hfiHB);   //  compute reference value from LISA + WD
    }
    
    if ((wdFLAG == 2) || (wdFLAG == 3))
    {
        //  look for frequency of interest in wd curve
        if ((fgw > fHB[0]) &&  (fgw < fHB[nHBCURVE-1]) )
        {
            jj = 0;
            while (fHB[jj] <= fgw) jj++;
            
            //  straight line fit on log-log graph
            
            m = (log10(hfHB[jj]/hfHB[jj-1]))/(log10(fHB[jj]/fHB[jj-1])); //  Slope btwn bkgnd pts
            b = log10(hfHB[jj]) - m*log10(fHB[jj]);    //  Intercept for fit btwn bkgnd pts
            
            hfiHB = pow(10.0,m*log10(fgw) + b);         //  Interpolated LISA value
            
            hfi = sqrt(hfiL*hfiL + hfiHB*hfiHB);   //  compute reference value from LISA + WD
        }
    }
    
    return hfi;   // return estimated value of LISA sensitivity curve at requested frequency
}

// =============================================================================
// =============================================================================


// ***************************************************************************
//
//  FUNCTION:  parseTokenInt
//
// ***************************************************************************

int parseTokenInt(char *myRead, char *delimeter)
{
    /* =============== VARIABLES =============== */
    
    char *tokenParse, *parseThis;  // strings for parsing tokens (elements) from a CSV line
    
    int inputNumber;
    
    /* =============== START ROUTINE HERE =============== */
    
    parseThis = myRead;                        // set up string for token parsing
    
    tokenParse = strsep(&parseThis,delimeter);    // parse up to the delimeter
    tokenParse = strsep(&parseThis,delimeter);    // parse after the delimeter
    
    sscanf(tokenParse,"%d", &inputNumber);       // return the string that comes after the DELIMETER
    
    return inputNumber;                           // return number from header line
    
} // FUNCTION: headerTokenInt


// ***************************************************************************
//
//  FUNCTION:  parseTokenLong
//
// ***************************************************************************

long parseTokenLong(char *myRead, char *delimeter)
{
    /* =============== VARIABLES =============== */
    
    char *tokenParse, *parseThis;  // strings for parsing tokens (elements) from a CSV line
    
    long inputNumber;
    
    /* =============== START ROUTINE HERE =============== */
    
    parseThis = myRead;                        // set up string for token parsing
    
    tokenParse = strsep(&parseThis,delimeter);    // parse up to the delimeter
    tokenParse = strsep(&parseThis,delimeter);    // parse after the delimeter
    
    sscanf(tokenParse,"%ld", &inputNumber);       // return the string that comes after the DELIMETER
    
    return inputNumber;                           // return number from header line
    
} // FUNCTION: headerTokenLong


// ***************************************************************************
//
//  FUNCTION:  parseTokenDouble
//
// ***************************************************************************

double parseTokenDouble(char *myRead, char *delimeter)
{
    /* =============== VARIABLES =============== */
    
    char tmpString[MAX_LENGTH + 1];
    char *tokenParse, *parseThis;  // strings for parsing tokens (elements) from a CSV line
    
    double inputNumber;
    
    /* =============== START ROUTINE HERE =============== */
    
    parseThis = myRead;                        // set up string for token parsing
    
    tokenParse = strsep(&parseThis,delimeter);    // parse up to the delimeter
    tokenParse = strsep(&parseThis,delimeter);    // parse after the delimeter
    
    sscanf(tokenParse,"%lf", &inputNumber);       // return the string that comes after the DELIMETER
    
    return inputNumber;                           // return number from header line
    
} // FUNCTION: headerTokenDouble



// ***************************************************************************
//
//  FUNCTION:  parseTokenString
//
// ***************************************************************************
void parseTokenString(char *myString, char *myRead, char *delimeter)
{
    /* =============== VARIABLES =============== */
    
    char *tokenParse, *parseThis;  // strings for parsing tokens (elements) from a CSV line
    
    /* =============== START ROUTINE HERE =============== */
    
    parseThis = myRead;                           // set up string for token parsing
    
    tokenParse = strsep(&parseThis,delimeter);    // parse up to the delimeter
    tokenParse = strsep(&parseThis,delimeter);    // parse after the delimeter
    
    sscanf(tokenParse,"%s", myString);            // return the string that comes after the DELIMETER
    
} // FUNCTION: headerTokenString


// =============================================================================
// =============================================================================


// ***************************************************************************
//
//  FUNCTION:  ErrorExit
//
//  General routine for error report and exit of program. Takes name of routine,
//  error message, and reports to stdout with the time of the exception, then
//  exits.
//
// ***************************************************************************

void ErrorExit(char routine[], char errorMsg[])
{
    /* =============== VARIABLES =============== */
    
    // --- date & time processing ---
    time_t errorTime;
    struct tm *errorQuery;
    char timeEnded[128];
    

    /* =============== START ROUTINE HERE =============== */

    // --- get the end time of the run ---
    time(&errorTime);                                 // get time from the system
    errorQuery = localtime(&errorTime);               // process time
    strftime(timeEnded,128,"%x - %I:%M%p",errorQuery);  // human readable string
    

    printf("\n\n");
    printf("============== ERROR & EXIT ===============\n\n");
    printf("ROUTINE   : %s\n",routine);
    printf("PROBLEM   : %s\n\n",errorMsg);
    printf("EXIT TIME : %s\n\n",timeEnded);
    printf("===========================================\n\n");

    fflush(stdin);
    exit(0);
    
    return;
    
} // FUNCTION: Error Exit
