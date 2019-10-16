// NEW MULTIPECHO SEQUENCE WITH PACE, WIP PACKAGE


//-----------------------------------------------------------------------------
//  Copyright (C) Siemens AG 1999  All Rights Reserved.  Confidential
//-----------------------------------------------------------------------------
//
// Project: NUMARIS/4
//
//    File: \n4_servers1\pkg\MrServers\MrImaging\seq\a_ep2d.cpp
//
//  Author: Clinical
//          Thomas Kluge;    Siemens AG Med MRIA/Seq;  (09131) 84-8049 (01/1999-02/2002)
//          Michael Zwanger; Siemens AG Med MREA-fMRI; (09131) 84-2672 (only ep2d_diff)
//          Josef Pfeuffer;  Siemens AG Med MR PLM AW Neuro; (only ep2d_pasl)
//
//    Lang: C++
//
// Descrip: Source Code for Standard NUMARIS/4 Sequences of type EPI.
//
// EGA Requirement Key: As shown on the following lines:
//
//   Abbrev.   Translation                                        Relevant for
//   -------   -----------                                        ------------
//   EGA-All   One, some or all of the following keys:            Any EGA requirement
//   EGA-01    {:IMPLEMENT:000_EGA_BildOri_SW_SequenzROVz::}      GR/GP   polarity
//   EGA-02    {:IMPLEMENT:000_EGA_BildPos_SW_SequenzSSelVz::}    GS      polarity
//   EGA-03    {:IMPLEMENT:000_EGA_BildMass_SW_SequenzROPC::}     GR/GP   amplitude
//   EGA-04    {:IMPLEMENT:000_EGA_BildPos_SW_SequenzSSel::}      GS      amplitude
//   EGA-05    {:IMPLEMENT:000_EGA_BildPos_SW_NCOFrequenzSSel::}  SRF     frequency
//   EGA-06    {:IMPLEMENT:000_EGA_BildPos_SW_NCOFrequenzRO::}    Readout frequency
//   EGA-07    {:IMPLEMENT:000_EGA_BildOri_SW_OrientierungTest::} Image orientation
//   EGA-11    {:IMPLEMENT:022_EGA_AbbFehler_SW_SpinEchoSequ.:}   GS ratio of refocussing and excitation pulses in 2D SE
//
//-----------------------------------------------------------------------------

/** \requirement N4_elh_DIFF_FLAIR,
	N4_elh_DTI_measurement_Parameter_Range_Matrix_PAT 
*/

//#define SEQSIGNATUREID "3000489643"     // ID needed for pTX !

//-------------------------------------------------------------------------------------
// imports
//-------------------------------------------------------------------------------------
#include      "MrServers/MrMeasSrv/SeqIF/sde_allincludes.h"
#include      "MrServers/MrImaging/libSeqUtil/ReorderInfoEPI.h"
#include      "MrServers/MrImaging/libSBB/SBBBinomialPulses.h"
// benpos WIP
#include      "SBBEPIKernel_WIP.h"

#include      "MrServers/MrImaging/seq/epi_StdUILink.h"
#include      "MrServers/MrMeasSrv/MeasUtils/MeasMath.h"  // minimum/maximum
#include      "MrServers/MrImaging/seq/a_ep2d.h"
#include      "MrServers/MrImaging/seq/epi_Config.h"
#include      "MrServers/MrPerAns/PerProxies/GCProxy.h"         // .getType
#include      "MrServers/MrImaging/seq/common/SeqLoopLongTRTrig/SeqLoopLongTRTrig.h"


#ifdef SUPPORT_iPAT_a_ep2d
//  #pragma message ("NOTE: This single shot EPI sequence supports iPAT.")
    #include  "MrServers/MrImaging/seq/common/iPAT/iPAT.h"
#endif

#ifndef VXWORKS
    #include  "MrServers/MrProtSrv/MrProtocol/libUILink/UILinkArray.h"
    #include  <vector>
#endif
	
#ifdef EP2D_DIFF
    #include "SBBDiffusion.h"
    #include "MrServers/MrImaging/libINIAccess/INIAccess.h"
#endif

#ifdef PASL  
    #include "SBBPasl.h"
	#include "MrServers/MrProtSrv/MrProt/Asl/Asl.h" //benpos VASO
#endif

#ifdef PACE3D
	// benpos WIP - include these files from the local seq directory
	// they havent been changed (aprt from the include paths)
	// but it seems very silly to incude them from \seq\a_ep2d_pace 
    #include      "PaceFeedback.h"
    #include      "PACE_fSEQReceive.h"
    
    // the global instance of the PACE feedback class
    extern PaceFeedback gPaceFeedback;
#endif



// benpos WIP 
	// special card include and globals 
	#include "parameter_map3beta.h"
	long	l_LongEntry_delayBetweenFirstAndSecondEcho =0;
	long	l_LongEntry_delayBetweenRemainingEchoes = 0;
	long    l_LongEntry_FixedZshim = 0; 
	long    l_LongEntry_AlternZshim = 0; 
	bool	b_Checkbox_LogPhysio = false;
	selection l_SelectionBox_EchoCombineMode;
	long	l_LongEntry_T2starweight = 30;
	bool    b_Checkbox_SaveSeparateEchoes = false;

	selection l_SelectionBox_PhaseCorrMode;

	selection l_SelectionBox_RefScanMode = 0;
	#include "MrServers/MrImaging/libSBB/SBBPATRefScan.h"
	static SeqBuildBlockPATRefScan SBBPATRefScan(NULL);
	bool b_sbbpatrefscans_done = false;

	long  l_LongEntry_FleetDummies = 5;
	long  l_LongEntry_FleetFlipAngle = 20;
	

//	double  d_DoubleEntry_SSBPatRefscanTE = 4;
//	double  d_DoubleEntry_SSBPatRefscanBW = 300;
//	double  d_DoubleEntry_SSBPatRefscanBaseRes = 64;

	long l_LongEntry_FatSatFlipAngle  = 110;
	
	// ggg changed to long double l_LongEntry_RfDuration = 2560;
	long l_LongEntry_RfDuration = 2560;
	double d_DoubleEntry_RfBWTP = 5.2;
	double d_DoubleEntry_FFTscale = 1; 
	
	// benpos SMS-CAIPI -->
	long l_LongEntry_SmsFactor = 1;
	long l_LongEntry_SmsShiftFactor = 1; 
	bool b_Checkbox_SmsPhaseOptimise = true; 
	bool b_Checkbox_SmsOnlineRecon =  true;
	// <-- benpos SMS-CAIPI
	
	// for writing protocol to file
	#ifndef VXWORKS
		char *env_custom_seq = getenv("CustomerSeq");
		char protfilename [256];
	#endif	
	
//benpos Pasl parameters
	long  l_LongEntry_trfociAmpl  = 100.0; // Renzo
	long  l_LongEntry_trfociBWDTH = 100.0; // Renzo
	long  l_LongEntry_trfociD = 100.0; // Renzo


	//benpos VASO parameter
	long l_VasoBoldFilltime = 0;



	long	l_LongEntry_VolumesPerInversion = 1;


	// for physio logging
	// use CPmuSequenceGuard class developed by Joost Kuijer at the Free University Amsterdam Medical Centre
	// this ensures that logging is ALWAYS stopped also in case
	// of sequence crash
	#include "MrServers\MrMeasSrv\PMU\pmusequence.h"
	#include "VUMC_CPmuSequenceGuard.h"
	#ifdef VXWORKS
		static CPmuSequenceGuard MyPhysioPuls;
		static CPmuSequenceGuard MyPhysioResp;
		static CPmuSequenceGuard MyPhysioECG;
		static CPmuSequenceGuard MyPhysioEXT;
	#endif


	//for better ghost correction in long scans, adapt the PC-algorythm (re. discussions with Thorsten Feiweier)
	#ifndef ICE_ONLINEPC_AUTOCROSSCORR_ACROSSSEGMENTS
		#define ICE_ONLINEPC_AUTOCROSSCORR_ACROSSSEGMENTS 64
	#endif
	
// ggg
#include "MrServers/MrImaging/seq/Spiral_Infra/gCommon.h"
#include "MrServers/MrImaging/seq/Spiral_Infra/gSpiral.h"

// parameter_map
static long paramLongROSamples = 7128; 			  // RO Samples
/*static long paramLongInterleaves = 1;				// Interleaves
static long paramLongSpiralType = 0;			   // Spiral Type
static long paramLongSpSlewRate=125; // for parameter map
static long paramLongSpGradAmp=25; // for parameter map
static long paramLongSpiralBW=400000; // for parameter map
static long paramLongZCToFlat = 0;*/

static double paramLongInterleaves = 1;				// Interleaves
static double paramLongSpiralType = 0;			   // Spiral Type
static double paramLongSpSlewRate=155; // for parameter map
static double paramLongSpGradAmp=35; // for parameter map
static double paramLongSpiralBW=400000; // for parameter map
static double paramLongZCToFlat = 0;

static double paramDoubleEchoShift = 0;
static double paramDoubleCAIPIPeriod = 0;
static double paramDoubleCAIPIDelay = 0;

static double paramDoubleVD = 0;
static double paramDoubleAccR = 1;

selection nIceProgram = 3; // Empty

#ifndef VXWORKS
/// Searches for the correct tool tip text for double parameters 
  unsigned _WIP_DOUBLE_GetToolTipVD(LINK_DOUBLE_TYPE* const pThis, char* arg_list[], long lIndex)
  {
    static char tToolTip[20000];
    
	sprintf(tToolTip,"%s",getSpiralToolTip(0));
	arg_list[0] = tToolTip;
	return MRI_STD_STRING;
  }

  unsigned _WIP_LONG_GetToolTipVD(LINK_LONG_TYPE* const pThis, char* arg_list[], long lIndex)
  {
    static char tToolTip[20000];
    
	sprintf(tToolTip,"%s",getSpiralToolTip(0));
	arg_list[0] = tToolTip;
	return MRI_STD_STRING;
  }
#endif // VXWORKS

//-------------------------------------------------------------------------------------
// sequence variant name
//-------------------------------------------------------------------------------------
#if   defined EP2D_DIFF
static const char *ptVariant = {"EPI 2D DIFFUSION Sequence"};
static const char *ptAddRTEB = {"EPI_Diffusion"};

#elif defined EP2D_FID
#ifdef PERF
static const char *ptVariant = {"N4 Standard EPI 2D FID Sequence (with perfusion postprocessing)"};
static const char *ptAddRTEB = {"nothing"}; 
#endif
#ifdef BOLD
// benpos WIP 
char *ptVariant = {"WIP 2D MULTIECHO EPI SEQUENCE Sequence (with no fmri postprocessing)"};
char *ptAddRTEB = {"nothing"};

#endif
#ifdef PASL
static const char *ptVariant = {"EPI 2D FID Sequence with Pulsed ASL preparation"};
static const char *ptAddRTEB = {"Pulsed ASL bipolar crusher gradient"};
#endif

#elif defined EP2D_SE
static const char *ptVariant = {"N4 Standard EPI 2D SE Sequence"};
static const char *ptAddRTEB = {"refocusing RF-pulse with spoilers"};

#else
#error Sequence Type Undefined
#endif


//-------------------------------------------------------------------------------------
// macro for error-handling
//-------------------------------------------------------------------------------------
#define  mSBBErrGotoFinish(A, B) {lStatus = (A).getNLSStatus(); if (!pSeqLim->isContextPrepForBinarySearch()) TRACE_PUT2_NLS(TC_INFO, TF_SEQ, "%s: %s",ptModule, B,lStatus); goto FINISHED;}

//-------------------------------------------------------------------------------------
// Helper function for cooling pause calculation
//
// - EPIKernel has to be prepared
// - REOInfo has to be prepared
//-------------------------------------------------------------------------------------
long lCoolingPause (bool bDynamic = true);


// benpos caipi --> helper function for slice index calculation
long getSmsRecalculatedSliceIndex(MrProt* pMrProt, sSLICE_POS* pSLC, long lOrigSlcIndex, long lSmsFactor);
// <-- benpos caipi

//-------------------------------------------------------------------------------------
// prototype of local functions
//-------------------------------------------------------------------------------------
static NLS_STATUS fSEQRunKernel          (MrProt*,SeqLim*,SeqExpo*,long,long,long,long);
bool              calculateTRTIFillTimes (MrProt*, SeqLim*, SeqExpo*,long*,long*);

//-------------------------------------------------------------------------------------
// GPA-proxy (performance data ?, PE enabled ?)
//-------------------------------------------------------------------------------------
static GPAProxy theGPA;  

//-------------------------------------------------------------------------------------
// GC-proxy
//-------------------------------------------------------------------------------------
static GCProxy theGC ;

//-------------------------------------------------------------------------------------
// MSU-proxy (for field strength)
//-------------------------------------------------------------------------------------
static MSUProxy theMSU;



/*###########################################################################

			    C l a s s:   S e q l o o p E P 2 D

###########################################################################*/

//-------------------------------------------------------------------------------------
// CHARM 305893 : TR for single shot EPI sequences should always be physically correct!
//
// This requires new functionality of SeqLoop which is realized in a class derived from
// SeqLoop. The interface of SeqLoop was adapded so that we are able to realize the
// following things within this module:
//
// - TR check of libSeqUT is performed also between two measurements to guarantee that
//   the modifications made are also tested by the sequence UT.
// - The time for SUBFINI / SUBSTRT between two repetitions is taken into account for
//   TR and measuerement time calculations.
// - The time for a halt in the case of single-shot triggering is now correctly handled
//   within the TR-period concerning TR and measurement time calculations.
// - SeqLoop can be asked not to distribute the slices evenly over a TR-period but use
//   a certain TR-fill at the end of each concatenation. This time period can be
//   specified with a new protocol parameter. Choosing the  maximum value of this time
//   period for a given TR results in the measurement of all slices without any TR-fill
//   times between them.
//
// CHARM 354879 : Adaptive IR slice thickness for SE and DIFF variants
//
// If an inversion pulse is used with diffusion epi, the doubled inversion thickness
// can generate crosstalk issues (e.g. STIR: imperfect fat suppression) if SeqLoop
// dictates a nesting of inversions (e.g. invert slice 1 - invert slice 3 - scan
// slice 1 - scan slice 3). Thus, an inversion thickness identical to the slice 
// thickness seems reasonable. 
// However, an inversion pulse is also used for fluid attenuation (FLAIR). Due to
// inflow effects, an increased inversion thickness significantly enhances the quality
// of fluid attenuation. 
// As a compromise, a doubled inversion thickness will be only used if the inversion
// time exceeds 500ms, which is a reasonable indication of a FLAIR protocol.
//
// SeqLoopLongTRTrig:
// SeqLoopEP2D is derived from SeqLoopLongTRTrig. This class supports the possibility 
// of activating a 'Long TR Triggering' mode, in which the slices are organised into 
// diffrent slice groups, which are acquired during different RR intervals. The number 
// of slice groups (and hence the number of RR intervals per volume) is set using the 
// number of concatenations parameter. Long TR Triggering mode makes it possible to
// use a long TR and still acquire all slices during diastole. This can be useful 
// in diffusion imaging of the brain to avoid artefacts associated with CSF pulsation. 
// If the flag SUPPORT_PACE is defined, the class SeqLoopLongTRTrig is derived from 
// the class SLFB to enable support for navigator triggering. Note however that the 
// two techniques cannot be used at the same time. Long TR triggering mode is 
// activated for some sequence types in fSEQPrep(). 
// 
//-------------------------------------------------------------------------------------

class SeqLoopEP2D : public SeqLoopLongTRTrig
{
public:
    
    SeqLoopEP2D();
    virtual bool calcFillTimesOnly (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo, long lWantedTR=-1, long lTimeForOSC=-1);
    virtual void calcMeasurementTimeUsec (MrProt* pMrProt, SeqLim* pSeqLim);
    virtual double getTotalMeasTimeUsec (MrProt* pMrProt, SeqLim* pSeqLim);
    virtual void doClockInitTRBetweenRepetetitions  (MrProt *pMrProt, SeqLim *pSeqLim, SeqExpo *pSeqExpo, long lSliceIndex, long lEcho);
    virtual bool insertTRFillEnd(long lFillTime);
    virtual void setPutTriggerDelayIntoTR(bool bValue);
    virtual void setAdditionalTimeForExtraEBinTRUsec(long lAdditionalTimeForExtraEBinTRUsec);
    virtual long prepGetFSDuration(MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo);
    virtual long prepGetOptFSDuration(MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo);
	// benpos SMS-CAIPI -->
	virtual void setSlicesToMeasure(long val);
	virtual void setInnerSliceNumber(long val);
	virtual long getInnerSliceNumber(void);
	virtual long getTRFillEnd(void);
	// <-- benpos SMS-CAIPI
#ifdef SUPPORT_PACE
    virtual bool runFreeLoop(pSEQRunKernel pf, MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo, sSLICE_POS* pSlcPos, sREADOUT* psADC);
#endif    

	//benpos - allow variable FatSat flip angle 
	virtual void setCSatFlipAngle(double dFlipAngle);	

protected:

    bool m_bPutTriggerDelayIntoTR;
    long m_lAdditionalTimeForExtraEBinTRUsec;
};


SeqLoopEP2D::SeqLoopEP2D ()
    :m_bPutTriggerDelayIntoTR(false)
    ,m_lAdditionalTimeForExtraEBinTRUsec(0)
{
    static const char * const ptModule = {"SeqLoopEP2D"};
    if (0) cout << ptModule; // to avoid warnings
}


void SeqLoopEP2D::setPutTriggerDelayIntoTR(bool bValue)
{
     m_bPutTriggerDelayIntoTR = bValue;
}


void SeqLoopEP2D::setAdditionalTimeForExtraEBinTRUsec(long lAdditionalTimeForExtraEBinTRUsec)
{
    m_lAdditionalTimeForExtraEBinTRUsec = lAdditionalTimeForExtraEBinTRUsec;
}


bool SeqLoopEP2D::calcFillTimesOnly (MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo,long lWantedTR_orig, long lTimeForOSC_orig)
{
    static const char * const ptModule = {"SeqLoopEP2D::calcFillTimesOnly"};

    //---------------------------------------------------------------------------
    // we expect that we are called from SeqLoop with the default parameters:
    //---------------------------------------------------------------------------
    if (lWantedTR_orig!=-1 || lTimeForOSC_orig!=-1)
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        textTR("unexpected ERROR: lWantedTR_orig!=-1 || lTimeForOSC_orig!=-1");
        setNLSStatus(SEQU_ERROR);
        return false;
    }

    //---------------------------------------------------------------------------
    // handling of m_lAdditionalTimeForExtraEBinTRUsec and m_TrigHaltDuration
    // may (or will) not work when multiple phases are selected
    //---------------------------------------------------------------------------
    if (m_PhasesToMeasure>1)
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        textTR("m_PhasesToMeasure>1 is currently not supported");
        setNLSStatus(SEQU_ERROR);
        return false;
    }

    //---------------------------------------------------------------------------
    // We want to force SeqLoop::calcFillTimesOnly to use a minimum for TRFillEnd
    // so we reduce the wanted TR for calculation of fill times.
    // 
    // For multiple measurements we want at least a TRFillEnd of 340us to be able
    // to compensate the increase of TR caused by the SUBFINI/SUBSTRT-event-block
    // inserted by SeqLoop. Compensation is done later by simply reducing the last
    // TRFilleEnd by 340us. Therefore TRFilleEnd must have at least this value. 
    //
    // Additional times must be put into TR:
    // - pMrProt->delayTimeInTR()
    // - m_lAdditionalTimeForExtraEBinTRUsec
    // - m_TrigHaltDuration
    //---------------------------------------------------------------------------
    long lWantedTR                     = pMrProt->tr()[0];
    long lReduceActualTRForCalculation = 0;

    lReduceActualTRForCalculation  = maximum(34L*GRAD_RASTER_TIME,pMrProt->delayTimeInTR());
    lReduceActualTRForCalculation += m_lAdditionalTimeForExtraEBinTRUsec;

    if (m_bPutTriggerDelayIntoTR && !m_TrigHaltSingleShot)
    {
        lReduceActualTRForCalculation += m_TrigHaltDuration;
    }

    lWantedTR -= lReduceActualTRForCalculation;

    if (lWantedTR<0)
    {
        lWantedTR=0;
    }


    //---------------------------------------------------------------------------
    // By increasing the time for the OSC-bit we force SeqLoop to put the ECG-fill
    // into TR, if m_TrigHaltSingleShot is true.
    //
    // Compare graphs of SeqLoop inner loop structure: fSBBECGFillTimeRun is called
    // directly before SBBOscBitRun.
    //---------------------------------------------------------------------------
    long lTimeForOSC = getScanTimeOscBit();

    if (m_bPutTriggerDelayIntoTR && m_TrigHaltSingleShot)
    {
        lTimeForOSC += getTrigHaltDuration();
    }

    //---------------------------------------------------------------------------
    // Do standard fill time calculation with modified wanted TR and time for OSC
    //---------------------------------------------------------------------------
    bool bRet = SeqLoopLongTRTrig::calcFillTimesOnly(pMrProt,pSeqLim,pSeqExpo,lWantedTR,lTimeForOSC);

    //---------------------------------------------------------------------------
    // correct members due to manipulations of lWantedTR
    //---------------------------------------------------------------------------
    for (long lI=0; lI<m_lConcatenations; ++lI)
    {
        m_SeqConcat[lI].m_TRFillEnd += lReduceActualTRForCalculation;
    }

    m_lActualTR += lReduceActualTRForCalculation;

    if (m_bHandleTRTIConflict && m_lTRneeded != 0)
    {
        m_lTRneeded += lReduceActualTRForCalculation;
    }

    //---------------------------------------------------------------------------
    // finished
    //---------------------------------------------------------------------------
    return bRet;
}


void SeqLoopEP2D::calcMeasurementTimeUsec (MrProt* pMrProt, SeqLim* pSeqLim)
{
    (void)SeqLoopLongTRTrig::calcMeasurementTimeUsec(pMrProt, pSeqLim);

    // we took care that ECG-fill-time was put into TR during calcFillTimesOnly,
    // standard SeqLoop-timing calculation does not recognize this, so we have to
    // do a correction here
    if (m_bPutTriggerDelayIntoTR)
    {
        m_dMeasureTimeInFirstMeasUsec  -= m_dTrigHaltTimeInFirstMeasUsec;
        m_dTrigHaltTimeInFirstMeasUsec  = 0;
    
        if (m_RepetitionsToMeasure)
        {
            m_dMeasureTimeInSecondMeasUsec  -= m_dTrigHaltTimeInSecondMeasUsec;
            m_dTrigHaltTimeInSecondMeasUsec  = 0;
        }
    }
}


double SeqLoopEP2D::getTotalMeasTimeUsec (MrProt* pMrProt, SeqLim* pSeqLim)
{
    (void)SeqLoopLongTRTrig::getTotalMeasTimeUsec(pMrProt, pSeqLim);

    // we take care that TR is always exact, so we consider the time for
    // SUBFINI / SUBSTRT between two repetitions during TR-calculation and during
    // execution of the sequence timing in calcFillTimesOnly
    m_dTotalMeasureTimeUsec -= 34.0 * GRAD_RASTER_TIME * m_RepetitionsToMeasure;

    #ifdef PASL  

      bool bEnableFirstPrepScanAsM0Scan = true;

      if( pMrProt->repetitions() <= 0 )
      {
        bEnableFirstPrepScanAsM0Scan = false;
      }
      
    #ifdef SUPPORT_iPAT_a_ep2d 
      
      // with iPAT a separate M0 scan BEFORE the iPAT scan is not possible
	    if( ( pMrProt->PAT().PATMode() != SEQ::PAT_MODE_NONE )
       && ( pMrProt->PAT().AccelFactPE() > 1               ) )
      {
        bEnableFirstPrepScanAsM0Scan = false;
      }
      
    #endif
      
      if( bEnableFirstPrepScanAsM0Scan )
      {
        m_dTotalMeasureTimeUsec += 34.0 * GRAD_RASTER_TIME;
      }

    #endif

    return m_dTotalMeasureTimeUsec;
}


void SeqLoopEP2D::doClockInitTRBetweenRepetetitions  (MrProt*, SeqLim*, SeqExpo*, long, long)
{
    // We never do a clock initilization for TR for the unit test between repetitions,
    // because TR must be physically correct for single shot EPI even over repetitions!
    return;
}


bool SeqLoopEP2D::insertTRFillEnd(long lFillTime)
{
    static const char * const ptModule = {"SeqLoopEP2D::insertTRFillEnd"};

    // SeqLoop will insert an event block for SUBFINI/SUBSTRT sync events which will disturb
    // our TR for multible measurements. So we have to shorten the last TR-fill-end of the
    // measurement to compensate this effect. 
    //
    if (     m_bIsLastScanInMeas
          && m_lRepetitionCounter<m_RepetitionsToMeasure
          && m_RepetitionsToMeasure
        )
    {
        lFillTime -= 340;
    }
    
    // adapt TRFillEnd for extra eventblocks in TR blocks
    if ( m_lRepetitionCounter<=m_RepetitionsToMeasure )
    {
      lFillTime -= m_lAdditionalTimeForExtraEBinTRUsec; 
    }

    // adapt TRFillEnd due to Trig-Halt
    if (m_bPutTriggerDelayIntoTR && !m_TrigHaltSingleShot)
    {
        lFillTime -= m_TrigHaltDuration; 
    }

    // check fill time:
    if (lFillTime<0)
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        textTR("FATAL: lFillTime<0");
        setNLSStatus(SEQU_ERROR);
        return false;
    }

    // setNLSStatus returns false for success
    return !setNLSStatus(fSBBFillTimeRun(lFillTime));
}

long SeqLoopEP2D::prepGetFSDuration(MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo)
{
    CSatFat.prep(pMrProt, pSeqLim, pSeqExpo);
    return CSatFat.getDurationPerRequest();
}

long SeqLoopEP2D::prepGetOptFSDuration(MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo)
{
    // The acutal dwell time of the OptFS module is calculated within 
    // calculateTRTIFillTimes, which requires a prepared SeqLoop. Since
    // we need the duration of the OptFS in advance, set dwell time to
    // zero here and calculate minimum OptFS duration (worst case with
    // respect to GPA load calculations). The correct dwell time will
    // be set at the end of fSeqPrep (when calculateTRTIFillTimes is
    // actually called).
    m_SBBOptfs.setDwellTime  (0);
    m_SBBOptfs.setJitterTime (0);
    m_SBBOptfs.prep(pMrProt, pSeqLim, pSeqExpo);
    return m_SBBOptfs.getDurationPerRequest();
}

#ifdef SUPPORT_PACE
bool SeqLoopEP2D::runFreeLoop (pSEQRunKernel pf, MrProt* pMrProt, SeqLim* pSeqLim, SeqExpo* pSeqExpo, sSLICE_POS* pSlcPos, sREADOUT* psADC)
{
    if( pMrProt->NavigatorParam().RespComp() == SEQ::RESP_COMP_OFF )
    {
        return SeqLoopLongTRTrig::runFreeLoop (pf,pMrProt,pSeqLim,pSeqExpo,pSlcPos,psADC);
    }
    //  The diffusion loop (direction and weighting) is inside the acquisition loop
    for (m_lOuterAcquisitionCounter = 0;m_lOuterAcquisitionCounter < m_AcquisitionsOuter; m_lOuterAcquisitionCounter++ ) 
    {
        psADC->Mdh.setCacq(static_cast<unsigned short>(m_lOuterAcquisitionCounter));
        // * -------------------------------------------------------------------------- *
        // * Free loop. Used e.g. for diffusion or Ciss                                 *
        // * -------------------------------------------------------------------------- *
        for (m_FreeLoopCounter = 0;m_FreeLoopCounter < m_FreeLoopLength; m_FreeLoopCounter++ ) 
        {
            if ( !runLineLoop(pf, pMrProt, pSeqLim, pSeqExpo, pSlcPos, psADC) )
            {
                TRACE_PUT2(TC_ALWAYS, TF_SEQ,"Error at %s(%d).",__FILE__,__LINE__);                
                return false;                
            }            
        } // End of free loop
    } // End Outer acquisition loop
    return true;  
}
#endif  //  SUPPORT_PACE

// benpos
void SeqLoopEP2D::setCSatFlipAngle(double dFlipAngle)
{
	// set FatSat flip angle
	CSatFat.setFlipAngle(dFlipAngle);
}

// benpos SMS-CAIPI -->
void SeqLoopEP2D::setSlicesToMeasure(long slices)
{
    m_SlicesToMeasure = slices;
}

void SeqLoopEP2D::setInnerSliceNumber(long innerslices)
{
    m_lInnerSliceNumber = innerslices;
}

long SeqLoopEP2D::getInnerSliceNumber(void)
{
    return m_lInnerSliceNumber;
}

long SeqLoopEP2D::getTRFillEnd(void)
{
    return m_SeqConcat[0].m_TRFillEnd;
}
// <-- benpos SMS-CAIPI


/*###########################################################################

			    G l o b a l   D e c l a r a t i o n s

###########################################################################*/

//-------------------------------------------------------------------------------------
// global variables
//-------------------------------------------------------------------------------------
static double      dMaxAllowedDataRate              = 1000.; // i.e. large
static bool        bPrepScansOnlyInFirstMeasurement = true;

#ifdef EP2D_DIFF
/// The diffusion loop counter runs over the b values and directions
/** It is required as a global variable because fSEQRunKernel() must export this
    information to fRunDIFFPlugInForEPIKernel(). The initialization with 0 is a dummy.
*/
static long               lDiffLoopCounter          = 0 ;
#endif

#ifdef SUPPORT_iPAT_a_ep2d
static long        lMinPrepScansNoPATRefScans       = 2;
// benpos SMS-CAIPI -->
//static long        lPrepScansNoPATRefScans          = 0;
static long        m_lPATPrepScans                  = 0;
static long        m_lInitialDummyScans              = 0;
static long        m_lPostDummyScans                = 0;
static long        m_lSliceAccPrepScans             = 0;
static long		   lPrepScanExcitationCounter = 0;  //FLEET
// <-- benpos SMS-CAIPI
static long        alPrepScanCounter[K_NO_SLI_MAX];
static bool        bSegmentedRefLines               = false;

#endif

//-------------------------------------------------------------------------------------
// standard loop structure realized within class
//-------------------------------------------------------------------------------------
static SeqLoopEP2D  mySeqLoop;

//-------------------------------------------------------------------------------------
// Slice position information (rotation matrices and shifts)
//-------------------------------------------------------------------------------------
static sSLICE_POS   asSLC[K_NO_SLI_MAX];

//-------------------------------------------------------------------------------------
// kernel calculation limits
//-------------------------------------------------------------------------------------
static KernelCalculationLimits myCalcLimits;

//-------------------------------------------------------------------------------------
// Sequence Building Blocks
//-------------------------------------------------------------------------------------
static SeqBuildBlockBinomialPulses   WE_121(NULL,false, 180.0, 3, "EPI_WE"); // the excitation SBB
static SeqBuildBlockEPIKernel        EPIKernel(NULL);                        // the EPI kernel class

// benpos SMS-CAIPI -->
static SeqBuildBlockExcitationRFPulse  SBBExcite(NULL);
static SeqBuildBlockExcitationRFPulse  SBBExciteFLEET(NULL);
static SeqBuildBlockExcitationRFPulse  SBBExciteSms(NULL);
// <-- benpos SMS-CAIPI


#ifdef PASL

// Fat sat all slices for Strong Fat Sat Mode
// We need to execute fat sat ourselves, instead of SeqLoop executing them because SeqLoop 
// will execute them before the PASL label pulses, and we want the fatsat after the label pulses
static SBBList                       SBB;                  // * List of SBBs         
static SeqBuildBlockCSat             CSatFat		( &SBB );  // * CSat fat suppression   
static SeqBuildBlockSpoilGrad        SpoilGrad	( &SBB );  // * Spoiler gradient after 
 
// group Pasl and RSat for display/positioning of labeling slab
static SBBList                       PaslSBB;
static SeqBuildBlockPasl             Pasl ( &PaslSBB );         
static SeqBuildBlockCrushGrad        PaslCrushGrad( NULL );
static bool                          bEnableFirstPrepScanAsM0Scan = false; //benpos true ->false
static long                          alPaslPrepScanCounter[K_NO_SLI_MAX];
static long                          alPaslScanCounter[K_NO_SLI_MAX];
#endif


//-------------------------------------------------------------------------------------
// reordering information data
//-------------------------------------------------------------------------------------
static ReorderInfoEPI REOInfo;

//-------------------------------------------------------------------------------------
// instance of EPI configuration class
//-------------------------------------------------------------------------------------
static epiConfig configInfo;

//-------------------------------------------------------------------------------------
// control execution of osc-bit
//-------------------------------------------------------------------------------------
#ifndef EP2D_DIFF
static const long lMaxOscBitSentFlags = 4096;
static bool abOscBitSentForMeas[lMaxOscBitSentFlags];
#endif

//-------------------------------------------------------------------------------------
// debug flags
//-------------------------------------------------------------------------------------
static long lDebug_SEQ_fSEQPrep = 0;


/*###########################################################################

				e p 2 d _ d i f f

###########################################################################*/

#ifdef EP2D_DIFF

/// This global object can be used to modify the sequence and SBB bahaviour for debugging
static INIAccess theConfig ("$CustomerSeq/ep2d_diff.ini") ;

/// This global SBBDiffusion object provides all diffusion functionality
static SBBDiffusion Diff ;



// ===========================================================================
//
//      fConfigureExcitationForEPIKernel
/*!
\brief  This function creates and configures the water excitation SBB used by the EPI Kernel

\return pointer to water excitation SBB

\author Michael.Zwanger@med.siemens.de
*/
// ===========================================================================

SeqBuildBlockExcitation* fConfigureExcitationForEPIKernel
(
    MrProt     *pMrProt, 
    SeqLim     *,
    SeqExpo    *pSeqExpo
)
{
    #undef  DEBUG_ORIGIN
    #define DEBUG_ORIGIN 0x00000200

    static const char *ptModule = {"fConfigureExcitationForEPIKernel"};

    static sRF_PULSE_EXT                 sSRF01("ExtExciteRF");

    sSRF01.setTypeExcitation       ();
    sSRF01.setFlipAngle            (pMrProt->flipAngle());
    sSRF01.setInitialPhase         (90.0);
    sSRF01.setFamilyName           ("SE2560A90.SE90_12A2_2");

    // The selection gradient amplitude of the excitation pulse has to be at least 20% different
    // from than that of the refocusing pulse in order to avoid ambiguity ('3rd arm') artefacts
    // (checked by SeqUT). Larger deviations help tp reduce signal contributions from off-resonance
    // spins (e.g. fat).
    // 
    // One has to make sure that the excitation pulses defined here match to the refocussing 
    // pulses defined in SBBDiffusion_Base.
    if ( theMSU.getNominalB0() < 2.0 )
    {
        sSRF01.setDuration             (2048);
    }
    else if ( theMSU.getNominalB0() < 4.0 )
    {
        sSRF01.setDuration             (2048);
    }
    else
    {
        sSRF01.setDuration             (4096);
    }

    sSRF01.setThickness            (pMrProt->sliceSeries().aFront().thickness());
    WE_121.setThickness            (pMrProt->sliceSeries().aFront().thickness());
   
    if(! WE_121.setExcitationRFPulse(&sSRF01,pMrProt,pSeqExpo))
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        textTR("WE_121.setExcitationRFPulse failed.");
        return NULL;
    }

    return &WE_121;
}



// ===========================================================================
//
//      fPrepDIFFPlugInForEPIKernel 
/*!
\brief  This function configures and prepares the SBBDiffusion

        This includes also to set the maximum gradient amplitude and slewrate.
        The RF pulses are prepared by preparing the SBBDiffusion.

        The gradient risetime must be set to a safe value which will not
        cause stimulations. If the GSWD detects possible stimulations, only
        the slew rates of the EPI readout will be changed, but not not those
        of the diffusion SBB. So, there is no help possible if the GSWD bays
        due to stimulations caused by the diffusion SBB.

        This function is registered as a call-back function in fSEQInit and
        executed during fSEQPrep.

\param  *pMrProt                pointer to the sequence protocol
\param  *pSeqLim                pointer to the sequence limits buffer
\param  *pSeqExpo               pointer to the sequence exports buffer
\param  lEPIKernelTEContributionBeforeRTEBPlugIn 
                Input: The time between the center of the RF pulse and the
                beginning of the plug-in
\param  lEPIKernelTEContributionAfterRTEBPlugIn  
                Input: The time between the end of the plug-in and the
                echo position
\param  *plRTEBPlugInDurationPerRequest
                Export: The function returns the duration for one call of the 
                plug-in.\n
                \b Attention: This time might exceed the available TE or TR.
                It is to be taken as the time \e needed by this plug-in.
\param  *plRTEBPlugInTEContribution
                Export: The function returns the TE contribution for one call 
                of the plug-in.
\param  *pdEnergyPerRequest_Ws
                Export: The function returns the energy for one call of the 
                plug-in.
\param  *pRTEBPlugInInvertsMagnetization   
                Export: The function returns if the plug-in does invert the
                magnetization.

\return NLS_STATUS code

\author Michael.Zwanger@med.siemens.de
*/
// ===========================================================================

NLS_STATUS fPrepDIFFPlugInForEPIKernel
(
    MrProt  *pMrProt,
    SeqLim  *pSeqLim,
    SeqExpo *pSeqExpo,
    long     lEPIKernelTEContributionBeforeRTEBPlugIn,
    long     lEPIKernelTEContributionAfterRTEBPlugIn,
    long    *plRTEBPlugInDurationPerRequest,
    long    *plRTEBPlugInTEContribution,
    double  *pdEnergyPerRequest_Ws,
    bool    *pRTEBPlugInInvertsMagnetization    
)
{
    static const char *ptModule={"fPrepDIFFPlugInForEPIKernel"};

    // --------------------------------------------------------------------------
    // configure SBBDiffusion
    // --------------------------------------------------------------------------
    if (! Diff.create(pMrProt))     {
      TRACE_PUT1(TC_INFO, TF_SEQ, "%s: Diff.create failed!", ptModule);
      return SEQU_ERROR;
    }
    
    //    Diff->setIniAccess(&theConfig);
    theConfig.rehash() ;
    
    if ( theConfig.getIntDefault(ptModule, "LANDMARK", 0) )
      TRACE_PUT5 (TC_INFO, TF_SEQ, 
		  "%s: SpinPrepTimeus=%ld, ADCusTillEcho=%ld, NoiseThreshold=%ld, expected SBB_plugin_duration=%ld",
		  ptModule, 
                  lEPIKernelTEContributionBeforeRTEBPlugIn, 
		  lEPIKernelTEContributionAfterRTEBPlugIn, 
		  pMrProt->diffusion().noiseLevel(), 
		  (pMrProt->te()[0] - lEPIKernelTEContributionAfterRTEBPlugIn - lEPIKernelTEContributionBeforeRTEBPlugIn) );
    
    Diff->setSpinPrepTimeus(lEPIKernelTEContributionBeforeRTEBPlugIn);
    Diff->setADCusTillEcho (lEPIKernelTEContributionAfterRTEBPlugIn);
    Diff->setNoiseThreshold(pMrProt->diffusion().noiseLevel());

    // Diffusion module needs some information on the GPA load of the readout module - if this information
    // is not provided, a default diffusion gradient amplitude similar to previous software versions is used.
    //
    // IMPORTANT: Do not rely on anything that is supposed to be calculated within SeqLoop here - SeqLoop
    // is not yet prepared!!
    {
        GPABalance theROBalance (&theConfig);

        // Check whether GPA is supported - otherwise, don't provide GPALoad and use 'old style' diffusion amplitudes
        if ( !theROBalance.lGetStatus() )
        {
            long int lI           = 0;
            long int lTimeUS      = 0;
            double   dSign        = 1.;
            double   dROAmp       = EPIKernel.getGRO().getAmplitude();      // We are lucky: the readout is prepared before this plugin
            long     lRORampUp    = EPIKernel.getGRO().getRampUpTime();
            long     lRORampDown  = EPIKernel.getGRO().getRampDownTime();
            long     lROFlatTop   = EPIKernel.getGRO().getFlatTopTime();
            long     lEchoSpacing = EPIKernel.getEchoSpacing();
            long     lEchoNumber  = REOInfo.getEchoTrainLength();

            if ( theConfig.getIntDefault(ptModule, "LANDMARK", 0) )
            {
                TRACE_PUT7 (TC_INFO, TF_SEQ, "%s: RO gradient parameters: Amp %7fmT/m, RampUp %ldus, FlatTop %ldus, RampDown %ldus, EchoSpacing %ldus, EchoTrain %ld", 
                    ptModule, dROAmp, lRORampUp, lROFlatTop, lRORampDown, lEchoSpacing, lEchoNumber);
            }
            
            // Add each readout gradient (with alternating sign)
            for (lI = 0; lI < lEchoNumber; ++lI)
            {
                theROBalance.lAddGradient (lTimeUS, dSign * dROAmp, lRORampUp, lRORampDown, lROFlatTop, GPABALANCE_X_AXIS);
                
                // Alternate the sign of the readout gradient
                dSign *= -1.;

                // Time increment
                lTimeUS += lEchoSpacing;
            }

            // Add delay: take into accout the additional cooling pause
            // Cooling pause is not considered here if IR or multiple concats (PACE) are enabled.
            // Background (see SeqLoopEP2D::calcFillTimesOnly above): in those cases, the whole
            // cooling pause is applied only after all kernel calls of the outer loop have taken
            // place. Thus, several kernels will in general be applied successively without any
            // forced cooling between them - still they should run!
            if ( (pMrProt->preparationPulses().inversion() == SEQ::INVERSION_OFF) && (pMrProt->concatenations() <= 1) )
            {
                // Do NOT consider the dynamic increase of the cooling pause here - the dynamic
                // increase is a RESULT of the GPA load calculation!
                //
                // Strategy: invest the reduction of the GPA load due to the STATIC (argument = false) 
                // pause into higher diffusion gradient amplitudes / shorter TE. The diffusion SBB then
                // calculates an additional DYNAMIC pause that is required for running the protocol
                // infinitely.
                lTimeUS += lCoolingPause (false);
            }

            // Add delay: RF excitation
            lTimeUS += WE_121.getDurationPerRequest();

            // Add delay: fat saturation
            //
            // We need prepared CSat's in order to get their time.
            // Quick and dirty: prepare them here - they will get
            // prepared again in mySeqLoop.prep() below ...
            switch ( pMrProt->preparationPulses().fatSuppression() )
            {
            case SEQ::FAT_SUPPRESSION_OPTIMAL:
                lTimeUS += mySeqLoop.prepGetOptFSDuration(pMrProt, pSeqLim, pSeqExpo);
                break;
            case SEQ::FAT_SATURATION:
                lTimeUS += mySeqLoop.prepGetFSDuration(pMrProt, pSeqLim, pSeqExpo);
                break;
            default:
                break;
            }

            // Add three phase correction echoes
            for (lI = 0; lI < 3; ++lI)
            {
                theROBalance.lAddGradient (lTimeUS, dSign * dROAmp, lRORampUp, lRORampDown, lROFlatTop, GPABALANCE_X_AXIS);
                
                // Alternate the sign of the readout gradient
                dSign *= -1.;

                // Time increment
                lTimeUS += lEchoSpacing;
            }

            // TRACE_PUT3 (TC_INFO, TF_SEQ, "%s: Cooling pause %li, RF duration %li", ptModule, lCoolingPause(false), WE_121.getDurationPerRequest() );

            // Provide diffusion module with the readout gradient events
            Diff->setGPALoad ( theROBalance );
        }
    }

    
    // --------------------------------------------------------------------------
    // prepare SBBDiffusion
    // --------------------------------------------------------------------------
    if (! Diff->prep (pMrProt, pSeqLim, pSeqExpo)) {
      if (!pSeqLim->isContextPrepForBinarySearch()) 
        TRACE_PUT1_NLS(TC_INFO, TF_SEQ, "%s: Diff->prep failed with ", ptModule, Diff->getNLSStatus());
      return Diff->getNLSStatus();
    }


    // --------------------------------------------------------------------------
    // set exports
    // --------------------------------------------------------------------------
    if ( plRTEBPlugInDurationPerRequest==NULL || plRTEBPlugInTEContribution     ==NULL ||
         pdEnergyPerRequest_Ws         ==NULL || pRTEBPlugInInvertsMagnetization==NULL  )   {
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s ERROR: NULL pointer to export parameters", ptModule);
        return SEQU_ERROR;
    }

    *plRTEBPlugInDurationPerRequest  = Diff->getDurationPerRequest();
    *plRTEBPlugInTEContribution      = Diff->getDurationPerRequest();
    *pdEnergyPerRequest_Ws           = Diff->getEnergyPerRequest();
    *pRTEBPlugInInvertsMagnetization = Diff->isMagnetizationInverted() ; 

    return SEQU_NORMAL;
}


// ===========================================================================
// 
//      fRunDIFFPlugInForEPIKernel 
/*!
\brief  This function runs SBBDiffusion::run

        This function is registered as a call-back function in \ref fSEQInit and
        executed during \ref fSEQRun.

\param  *pMrProt                pointer to the sequence protocol
\param  *pSeqLim                pointer to the sequence limits buffer
\param  *pSeqExpo               pointer to the sequence exports buffer
\param  *pSLC                   pointer to the actual slice
\param  lRTEBPlugInTEFillBefore is ignored. Should be zero.
\param  lRTEBPlugInTEFillAfter  is ignored. Should be zero.

\return NLS_STATUS code

\author Michael.Zwanger@med.siemens.de
*/
// ===========================================================================

NLS_STATUS fRunDIFFPlugInForEPIKernel
(
    MrProt      *pMrProt,
    SeqLim      *pSeqLim,
    SeqExpo     *pSeqExpo,
    sSLICE_POS  *pSLC,
    long         lRTEBPlugInTEFillBefore,
    long         lRTEBPlugInTEFillAfter
)
{
    static const char *ptModule={"fRunDIFFPlugInForEPIKernel"};


    // --------------------------------------------------------------------------
    // SBBDiffusion completely consumes the TE-fill times available, also it has
    // no methods to specify additional fill times.
    // So, if the RTEBPlugIn is told to insert fill times, we have a problem.
    // --------------------------------------------------------------------------
    if (lRTEBPlugInTEFillBefore!=0)
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s WARNING: lRTEBPlugInTEFillBefore=%ld (should be 0).", 
		   ptModule, lRTEBPlugInTEFillBefore) ;
    if (lRTEBPlugInTEFillAfter!=0)
        TRACE_PUT2(TC_INFO, TF_SEQ, "%s WARNING: lRTEBPlugInTEFillAfter=%ld (should be 0).", 
		   ptModule, lRTEBPlugInTEFillAfter) ;

    if (! Diff->run(pMrProt, pSeqLim, pSeqExpo, pSLC, lDiffLoopCounter, EPIKernel.getReadOutAddress() ) )  {
      // lDiffLoopCounter is a global variable which has been set in fSeqRunKernel
      TRACE_PUT1_NLS(TC_INFO, TF_SEQ, "%s: Diff->run failed with ", ptModule, Diff->getNLSStatus());
      return Diff->getNLSStatus();
    }
    
    return SEQU_NORMAL;
}
#endif





/*###########################################################################

				e p 2 d _ f i d

###########################################################################*/

#ifdef EP2D_FID
/*[ Function ****************************************************************\
*
* Name        : fConfigureExcitationForEPIKernel
*               
* Description : Creates and configures the SBB for excitation used by the EPI
*               Kernel.
*
* Return      : pointer to excitation SBB
*
\****************************************************************************/
SeqBuildBlockExcitation* fConfigureExcitationForEPIKernel
(
    MrProt     *pMrProt, 
    SeqLim     *,
    SeqExpo    *pSeqExpo
)
{
    #undef  DEBUG_ORIGIN
    #define DEBUG_ORIGIN 0x00000200

    static const char *ptModule = {"fConfigureExcitationForEPIKernel"};
/*
    static sRF_PULSE_SINC                sSRF01("SincRFPulse");

    sSRF01.setTypeExcitation      ();
    sSRF01.setFlipAngle           (pMrProt->flipAngle());
    sSRF01.setInitialPhase        (90);
    sSRF01.setBandwidthTimeProduct(5.2);
    sSRF01.setDuration            (2560);
    sSRF01.setSamples             ( 512);


	//benpos WIP
	sSRF01.setDuration            (l_LongEntry_RfDuration);
	sSRF01.setSamples             (l_LongEntry_RfDuration/5);
	sSRF01.setBandwidthTimeProduct(d_DoubleEntry_RfBWTP);

    sSRF01.setThickness           (pMrProt->sliceSeries().aFront().thickness());
    WE_121.setThickness           (pMrProt->sliceSeries().aFront().thickness());
   
    if(! WE_121.setExcitationRFPulse(&sSRF01,pMrProt,pSeqExpo))
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        textTR("WE_121.setExcitationRFPulse failed.");
        return NULL;
    }

    return &WE_121;
	
	
	*/	
	
	// benpos SMS-CAIPI -->
	
	//all the relevant pulse calculations happen in here in the include file:
	#include "make_multiplex_RF_ARB_WIP.h"

	if(! SBBExcite.setExcitationRFPulse(&sSRF_TEMP,pMrProt,pSeqExpo))
	{
		textTR("SBBExcite.setExcitationRFPulse failed.");
		return NULL;
	}
	if( !SBBExciteSms.setExcitationRFPulse(&sSRF_MPLX,pMrProt,pSeqExpo))
	{
		textTR("SBBExcite.setExcitationRFPulse failed.");
		return NULL;
	}

	//we feed SMS excitation to the kernel ourselves...
	EPIKernel.setSBBExciteSms(&SBBExciteSms); 

    // make a copy of the single-band pulse for FLEET
//	static sRF_PULSE_EXT    sSRF_FLEET("ExtExcFleet");
//	sSRF_FLEET = sSRF_TEMP;
//    sSRF_FLEET.setFlipAngle            (double ( l_LongEntry_FleetFlipAngle ) );
static sRF_PULSE_SINC                sSRF01("SincRFPulse");

    sSRF01.setTypeExcitation      ();
    sSRF01.setFlipAngle           (double (l_LongEntry_FleetFlipAngle) );
    sSRF01.setInitialPhase        (90);

	//benpos WIP
	sSRF01.setDuration            (l_LongEntry_RfDuration);
	sSRF01.setSamples             (l_LongEntry_RfDuration/5);
	sSRF01.setBandwidthTimeProduct(5.2);

    sSRF01.setThickness           (pMrProt->sliceSeries().aFront().thickness());
    
	
	if( !SBBExciteFLEET.setExcitationRFPulse(&sSRF01,pMrProt,pSeqExpo))
	{
		textTR("SBBExciteFLEET.setExcitationRFPulse failed.");
		return NULL;
	}
	// .. end feed it to the EPI kernel.
	EPIKernel.setSBBExciteFLEET(&SBBExciteFLEET); 

	return &SBBExcite;
	// <-- benpos SMS-CAIPI
	
	
}
#endif





/*###########################################################################

        P A S L

###########################################################################*/

#ifdef PASL

/*[ Function ****************************************************************\
*
* Name        : fPrepPaslPlugInForEPIKernel
*               
* Description : prepares bipolar crusher gradient
*               
* Return      : NLS_STATUS code
*
\****************************************************************************/
NLS_STATUS fPrepPaslPlugInForEPIKernel
(
    MrProt  *pMrProt,
    SeqLim  *pSeqLim,
    SeqExpo *pSeqExpo,
    long     ,//lEPIKernelTEContributionBeforeRTEBPlugIn
    long     ,//lEPIKernelTEContributionAfterRTEBPlugIn
    long    *plRTEBPlugInDurationPerRequest,
    long    *plRTEBPlugInTEContribution,
    double  *pdEnergyPerRequest_Ws,
    bool    *pRTEBPlugInInvertsMagnetization    
)
{
    static const char *ptModule={"fPrepPaslPlugInForEPIKernel"};

    // --------------------------------------------------------------------------
    // prepare crusher
    // --------------------------------------------------------------------------
    PaslCrushGrad.setFlowLimit( Pasl.getFlowLimit(pMrProt) );
    if (! PaslCrushGrad.prep(pMrProt, pSeqLim, pSeqExpo))
    {
      TRACE_PUT1_NLS(TC_INFO, TF_SEQ, "%s: PaslCrushGrad.prep failed with ", ptModule, PaslCrushGrad.getNLSStatus());
      return PaslCrushGrad.getNLSStatus();
    }

    // --------------------------------------------------------------------------
    // set exports
    // --------------------------------------------------------------------------
    if (   plRTEBPlugInDurationPerRequest==NULL || plRTEBPlugInTEContribution     ==NULL
        || pdEnergyPerRequest_Ws         ==NULL || pRTEBPlugInInvertsMagnetization==NULL
       )
    {
        // caused by programming error
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s: at least one pointer to exports is NULL ",ptModule);
        return SEQU_ERROR;
    }

    *plRTEBPlugInDurationPerRequest   = PaslCrushGrad.getCrushTime(); 
    *plRTEBPlugInTEContribution       = PaslCrushGrad.getCrushTime();
    *pdEnergyPerRequest_Ws            = 0;
    *pRTEBPlugInInvertsMagnetization  = false;

    // --------------------------------------------------------------------------
    // if this point can be reached successfully
    // --------------------------------------------------------------------------
    return SEQU_NORMAL;
}

/*[ Function ****************************************************************\
*
* Name        : fRunPaslPlugInForEPIKernel
*               
* Description : executes bipolar crusher gradient
*               
* Return      : NLS_STATUS code
*
\****************************************************************************/
NLS_STATUS fRunPaslPlugInForEPIKernel
(
    MrProt      *pMrProt,
    SeqLim      *pSeqLim,
    SeqExpo     *pSeqExpo,
    sSLICE_POS  *pSLC,
    long         lRTEBPlugInTEFillBefore,
    long         lRTEBPlugInTEFillAfter
)
{
    static const char *ptModule={"fRunPaslPlugInForEPIKernel"};

    // --------------------------------------------------------------------------
    // run fill time before plugin
    // --------------------------------------------------------------------------
    if (lRTEBPlugInTEFillBefore) {
        fSBBFillTimeRun (lRTEBPlugInTEFillBefore);
    }

    // --------------------------------------------------------------------------
    // run crusher gradient
    // --------------------------------------------------------------------------
    if (! PaslCrushGrad.run(pMrProt,pSeqLim,pSeqExpo, pSLC))
    {
      TRACE_PUT1_NLS(TC_INFO, TF_SEQ, "%s: PaslCrushGrad.run failed with ", ptModule, PaslCrushGrad.getNLSStatus());
      return PaslCrushGrad.getNLSStatus();
    }

    // --------------------------------------------------------------------------
    // run fill time after plugin
    // --------------------------------------------------------------------------
    if (lRTEBPlugInTEFillAfter) {
      fSBBFillTimeRun (lRTEBPlugInTEFillAfter);
    }

    // --------------------------------------------------------------------------
    // if this point can be reached successfully
    // --------------------------------------------------------------------------
    return SEQU_NORMAL;
}

#endif  // PASL





/*###########################################################################

				e p 2 d _ s e

###########################################################################*/



#ifdef EP2D_SE

#include "MrServers/MrImaging/seq/epi_SERefocRTEB.h"
static   SERefocRTEB REFOCRFWithSpoilers;

/*[ Function ****************************************************************\
*
* Name        : fConfigureExcitationForEPIKernel
*               
* Description : Creates and configures the SBB for excitation used by the EPI
*               Kernel.
*
* Return      : pointer to excitation SBB
*
\****************************************************************************/
SeqBuildBlockExcitation* fConfigureExcitationForEPIKernel
(
    MrProt     *pMrProt, 
    SeqLim     *,
    SeqExpo    *pSeqExpo
)
{
    #undef  DEBUG_ORIGIN
    #define DEBUG_ORIGIN 0x00000200

    static const char *ptModule = {"fConfigureExcitationForEPIKernel"};

    static sRF_PULSE_EXT                 sSRF01("ExtExciteRF");

    sSRF01.setTypeExcitation       ();
    sSRF01.setFlipAngle            (pMrProt->flipAngle());
    sSRF01.setInitialPhase         (90.0);
    sSRF01.setFamilyName           ("SE2560A90.SE90_12A2_2");

    // Bandwidth of excitation pulse 20% higher than that of refocusing pulse to avoid 3rd arm artefact (checked by SeqUT)
    // Ensure that these values fit to those of the refocussing pulses defined in fPrepSEPlugInForEPIKernel
    if ( theMSU.getNominalB0() < 2.0 )
    {
        sSRF01.setDuration             (2048);
    }
    else if ( theMSU.getNominalB0() < 4.0 )
    {
        sSRF01.setDuration             (3072);
    }
    else
    {
        sSRF01.setDuration             (4096);
    }

    sSRF01.setThickness            (pMrProt->sliceSeries().aFront().thickness());
    WE_121.setThickness            (pMrProt->sliceSeries().aFront().thickness());
   
    if(! WE_121.setExcitationRFPulse(&sSRF01,pMrProt,pSeqExpo))
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        textTR("WE_121.setExcitationRFPulse failed.");
        return NULL;
    }

    return &WE_121;
}


/*[ Function ****************************************************************\
*
* Name        : fPrepSEPlugInForEPIKernel
*               
* Description : prepares refocussing RF pulse and its gradients
*               
* Return      : NLS_STATUS code
*
\****************************************************************************/
NLS_STATUS fPrepSEPlugInForEPIKernel
(
    MrProt  *pMrProt,
    SeqLim  *pSeqLim,
    SeqExpo *pSeqExpo,
    long     ,//lEPIKernelTEContributionBeforeRTEBPlugIn
    long     ,//lEPIKernelTEContributionAfterRTEBPlugIn
    long    *plRTEBPlugInDurationPerRequest,
    long    *plRTEBPlugInTEContribution,
    double  *pdEnergyPerRequest_Ws,
    bool    *pRTEBPlugInInvertsMagnetization    
)
{
    static const char *ptModule={"fPrepSEPlugInForEPIKernel"};

    // --------------------------------------------------------------------------
    // gradient performance depending on protocol
    // --------------------------------------------------------------------------
    REFOCRFWithSpoilers.m_dMaxMagnitude = theGPA.getGradMaxAmpl( SEQ::GRAD_FAST );

    if (pMrProt->gradSpec().isGSWDMode())
    {
        REFOCRFWithSpoilers.m_dMinRiseTime = pMrProt->gradSpec().GSWDMinRiseTime();
    }
    else
    {
        REFOCRFWithSpoilers.m_dMinRiseTime = theGPA.getGradMinRiseTimeAbsolute();
    }

    // --------------------------------------------------------------------------
    // configure refocusing RF
    // --------------------------------------------------------------------------
    static sRF_PULSE_EXT sSRF02("ExtRefocRF");

    sSRF02.setTypeRefocussing ();
    sSRF02.setFlipAngle       (180.0);
    sSRF02.setInitialPhase    (0.0);
    sSRF02.setThickness       (pMrProt->sliceSeries().aFront().thickness());
    sSRF02.setFamilyName      ("SE2560A180.SE180_12A2_2");
    if ( theMSU.getNominalB0() < 2.0 )
    {
    sSRF02.setDuration        (2560);
    }
    else if ( theMSU.getNominalB0() < 4.0 )
    {
        sSRF02.setDuration        (3840);
    }
    else
    {
        sSRF02.setDuration        (5120);
    }

    REFOCRFWithSpoilers.setRFPulseForRefocusing(&sSRF02);

    // --------------------------------------------------------------------------
    // configure calculation limits
    // --------------------------------------------------------------------------
    REFOCRFWithSpoilers.m_pCalcLimits = &myCalcLimits;

    // --------------------------------------------------------------------------
    // prepare REFOCRFWithSpoilers
    // --------------------------------------------------------------------------
    if (!REFOCRFWithSpoilers.prep(pMrProt,pSeqLim,pSeqExpo))
    {
        if (!pSeqLim->isContextPrepForBinarySearch()) TRACE_PUT1_NLS(TC_INFO, TF_SEQ, "%s: REFOCRFWithSpoilers.prep failed with ",ptModule, REFOCRFWithSpoilers.getNLSStatus());
        return REFOCRFWithSpoilers.getNLSStatus();
    }

    // --------------------------------------------------------------------------
    // set exports
    // --------------------------------------------------------------------------
    if (   plRTEBPlugInDurationPerRequest==NULL || plRTEBPlugInTEContribution     ==NULL
        || pdEnergyPerRequest_Ws         ==NULL || pRTEBPlugInInvertsMagnetization==NULL
       )
    {
        // caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s: at least one pointer to exports is NULL ",ptModule);
        return SEQU_ERROR;
    }

    *plRTEBPlugInDurationPerRequest  = REFOCRFWithSpoilers.m_lSBBDurationPerRequest_us;
    *plRTEBPlugInTEContribution      = REFOCRFWithSpoilers.m_lTEContribution;
    *pdEnergyPerRequest_Ws           = REFOCRFWithSpoilers.m_dEnergyPerRequest_Ws;
    *pRTEBPlugInInvertsMagnetization = true;    

    // --------------------------------------------------------------------------
    // if this point can be reached successfully
    // --------------------------------------------------------------------------
    return SEQU_NORMAL;
}

/*[ Function ****************************************************************\
*
* Name        : fRunSEPlugInForEPIKernel
*               
* Description : executes refocussing RF pulse and its gradients
*               
* Return      : NLS_STATUS code
*
\****************************************************************************/
NLS_STATUS fRunSEPlugInForEPIKernel
(
    MrProt      *pMrProt,
    SeqLim      *pSeqLim,
    SeqExpo     *pSeqExpo,
    sSLICE_POS  *pSLC,
    long         m_lRTEBPlugInTEFillBefore,
    long         m_lRTEBPlugInTEFillAfter
)
{
    static const char *ptModule={"fRunSEPlugInForEPIKernel"};

    // --------------------------------------------------------------------------
    // remaining configuration of REFOCRFWithSpoilers
    // --------------------------------------------------------------------------
    REFOCRFWithSpoilers.m_lStartTimeInEventBlock = m_lRTEBPlugInTEFillBefore;
    REFOCRFWithSpoilers.m_lFillEnd               = m_lRTEBPlugInTEFillAfter;

    // --------------------------------------------------------------------------
    // run REFOCRFWithSpoilers
    // --------------------------------------------------------------------------
    if (!REFOCRFWithSpoilers.run(pMrProt,pSeqLim,pSeqExpo, pSLC)) /*! EGA-07 !*/
    {
        TRACE_PUT1_NLS(TC_INFO, TF_SEQ, "%s: REFOCRFWithSpoilers.run failed with ",ptModule, REFOCRFWithSpoilers.getNLSStatus());
        return REFOCRFWithSpoilers.getNLSStatus();
    }

    // --------------------------------------------------------------------------
    // if this point can be reached successfully
    // --------------------------------------------------------------------------
    return SEQU_NORMAL;
}

#endif


/*###########################################################################

				e p 2 d 

###########################################################################*/

/*[ Function ****************************************************************\
*
* Name        : calculateTRTIFillTimes
*               
* Description : Calculates TR- and TI-fill times for SeqLoop and required TR
*               and TI values.
*
* Return      : true, if success
*
\****************************************************************************/
bool calculateTRTIFillTimes
(
    MrProt*  pMrProt,
    SeqLim*  pSeqLim,
    SeqExpo* pSeqExpo,
    long*    plNeededTI,
    long*    plNeededTR
)
{
    #undef  DEBUG_ORIGIN
    #define DEBUG_ORIGIN 0x00000200

    // ---------------------------------------------------------------------------
    // calculate the basic scan time and time needed for sats.  Needed for spir calculation
    // it is necessary to add the time spent for the spir pulse but first this must be calculated
    // --------------------------------------------------------------------------- 
    long lScanTimeSatsEtc = mySeqLoop.getlScanTimeAllSats();                      // SBB scan time
    long lScanTimeBasic   = EPIKernel.getDurationPerRequest() + lCoolingPause();  // Kernel scan time including mandatory fill time
    long lScanTimeStore   = lScanTimeSatsEtc;                                     // Helper



    // The sequence asks for a pointer to the SPAIR SBB from SeqLoop.  Using this pointer it then uses the SBB to calculate the required inversion time.
    // This is set in the SBB from within the SBB.
    if (pMrProt->preparationPulses().fatSuppression()==SEQ::FAT_SUPPRESSION_OPTIMAL )
    {
        SeqBuildBlockOptfs      * pSBBOptfs     = mySeqLoop.getpOptfs();
        SeqBuildBlockOptfsPrep  * pSBBOptfsPrep = mySeqLoop.getpOptfsPrep();

        // Add SPAIR time (depending on protocol parameters, e.g. TR) to lScanTimeSatsEtc
        if(pSBBOptfs == NULL || !pSBBOptfs->calcSPIRTime (pMrProt,pSeqLim,pSeqExpo,lScanTimeBasic,lScanTimeSatsEtc,0,SeqBuildBlockOptfs::SPIR_CALC_TYPE_EPI,pSBBOptfsPrep))
        {
            TRACE_PUT1(TC_ALWAYS, TF_SEQ,"%s: The calc SPIR time has failed, prob due to an error in the tickle pulse  ", ptModule);
            return false;
        }
    }

    // Update SBB scan time including SPAIR
    mySeqLoop.setlSBBScanTime (lScanTimeSatsEtc);
    
    #ifdef PASL


	//add fat sat time for strong mode
	if (pMrProt->preparationPulses().fatSatMode() == SEQ::FAT_SAT_STRONG)
	{
		//for regular ASL and VASO without BOLD
		if  ( pMrProt->Asl().Mode() != SEQ::ASL_CUSTOM2) 
		{
			mySeqLoop.setlKernelScanTime( (lScanTimeBasic + CSatFat.getDurationPerRequest() + SpoilGrad.getDurationPerRequest() ) * l_LongEntry_VolumesPerInversion ); //benpos
		}
		else
		{
			//for VASO with BOLD, add BOLD filltime / Numberof slices
			//VASO_TI1 --> pMrProt->ti()[0]  / 1000 ; 
			//BOLD_TI2 --> pMrProt->ti()[1]  / 1000 ; 
			while (l_VasoBoldFilltime < 0) 
			{ 
				pMrProt->ti()[1] = pMrProt->ti()[1] + 100;
				l_VasoBoldFilltime =  (pMrProt->ti()[1] - pMrProt->ti()[0] )  -  pMrProt->sliceSeries().size()/l_LongEntry_SmsFactor * ( lScanTimeBasic + CSatFat.getDurationPerRequest() + SpoilGrad.getDurationPerRequest() ) ;
			}
			l_VasoBoldFilltime =  (pMrProt->ti()[1] - pMrProt->ti()[0] )  -  pMrProt->sliceSeries().size()/l_LongEntry_SmsFactor * ( lScanTimeBasic + CSatFat.getDurationPerRequest() + SpoilGrad.getDurationPerRequest() ) ;
      
			// kernelscantime is for one slice, so here we add l_VasoBoldFilltime/Nslices
			mySeqLoop.setlKernelScanTime( (lScanTimeBasic + CSatFat.getDurationPerRequest() + SpoilGrad.getDurationPerRequest() ) * l_LongEntry_VolumesPerInversion  + l_VasoBoldFilltime/pMrProt->sliceSeries().size()*l_LongEntry_SmsFactor ); //benpos
cout << "l_VasoBoldFilltime " << l_VasoBoldFilltime << endl;		
		}

	}
	else
	{
		//for regular ASL and VASO without BOLD
		if  ( pMrProt->Asl().Mode() != SEQ::ASL_CUSTOM2) 
		{	
			mySeqLoop.setlKernelScanTime( lScanTimeBasic * l_LongEntry_VolumesPerInversion ); //benpos
		}
		else
		{
			//for VASO with BOLD, add BOLD filltime / Numberof slices
			//VASO_TI1 --> pMrProt->ti()[0]  / 1000 ; 
			//BOLD_TI2 --> pMrProt->ti()[1]  / 1000 ; 
			while (l_VasoBoldFilltime < 0) 
			{ 
				pMrProt->ti()[1] = pMrProt->ti()[1] + 100;
				l_VasoBoldFilltime =  (pMrProt->ti()[1] - pMrProt->ti()[0] )  -  pMrProt->sliceSeries().size()/l_LongEntry_SmsFactor * lScanTimeBasic  ;
			}
			l_VasoBoldFilltime =  (pMrProt->ti()[1] - pMrProt->ti()[0] )  -  pMrProt->sliceSeries().size()/l_LongEntry_SmsFactor * lScanTimeBasic  ;
						
			// kernelscantime is for one slice, so here we add l_VasoBoldFilltime/Nslices
			mySeqLoop.setlKernelScanTime( lScanTimeBasic  * l_LongEntry_VolumesPerInversion  + l_VasoBoldFilltime/pMrProt->sliceSeries().size()*l_LongEntry_SmsFactor ); //benpos
cout << "l_VasoBoldFilltime " << l_VasoBoldFilltime << endl;
		}
	}
    #else
		mySeqLoop.setlKernelScanTime(lScanTimeBasic);
    #endif

	// benpos SMS-CAIPI -->
	//for the following timing calculations, this tricks seqloop to calculate for slice accelerated TR
	mySeqLoop.setSlicesToMeasure( pMrProt->sliceSeries().size() / l_LongEntry_SmsFactor );
    // <-- benpos SMS-CAIPI
	

    bool bSuccess = mySeqLoop.calcFillTimes ( pMrProt, pSeqLim, pSeqExpo );


    if (!bSuccess) 
    {
        // What follows here is basically a copy of fUICEvaluateSeqLoopTrTiFillTimes. However, a special
        // handling is required in case that SPAIR is enabled: if the sequence demands to increase TR by
        // a certain amount, the SPAIR time needs to be updated which in turn might require an additional
        // TR increase. This hen-egg-problem needs to be explicitely resolved.

        if (mySeqLoop.getlTRneeded()!=0 || mySeqLoop.getlTIneeded()!=0)
        {
            // We could proceed, if we change TR and/or TI !
            
            if (mySeqLoop.getlTRneeded()!=0)
            {
                // NOTE: In contrast to TI mySeqLoop does NOT take care of the TR
                //       increment. So we have to do it here.

                double dInc = pSeqLim->getTR()[0].getInc();
                
                // Avoid division by zero
                if ( fabs(dInc) < 1.e-10 )
                {
                    // Unexpected error
                    // => trace also if pSeqLim->isContextPrepForBinarySearch()
                    TRACE_PUT0(TC_INFO, TF_SEQ, "ERROR in calculateTRTIFillTimes(): dInc == 0. ");
                    return false;
                }
                
                // ("..........NOTE!!!!! we assume something about TR-limit-handling in MrUILink !")
                if (mySeqLoop.getlTRneeded() >   10000) 
                {
                    dInc *= 10.0;
                }
                if (mySeqLoop.getlTRneeded() > 1000000) 
                {
                    dInc *= 10.0;
                }
                
                *plNeededTR = (long)(0.5+dInc*ceil(mySeqLoop.getlTRneeded()/dInc));
                
                
                // SPAIR: 
                // in this case we have to consider NeededTR as the actual TR (because it will be applied by the solve handler)
                // and thus perform the calculation again!
                if (pMrProt->preparationPulses().fatSuppression()==SEQ::FAT_SUPPRESSION_OPTIMAL )
                {
                    // Make a copy of MrProt and set the needed TR
                    MrProt sTmp(*pMrProt);
                    sTmp.tr()[0] = *plNeededTR;

                    // Revert original SBB scan time excluding SPAIR
                    lScanTimeSatsEtc = lScanTimeStore;
                    
                    SeqBuildBlockOptfs      * pSBBOptfs     = mySeqLoop.getpOptfs();
                    SeqBuildBlockOptfsPrep  * pSBBOptfsPrep = mySeqLoop.getpOptfsPrep();
                    // Add SPAIR time (depending on protocol parameters, e.g. TR) to lScanTimeSatsEtc
                    if(pSBBOptfs == NULL || !pSBBOptfs->calcSPIRTime (&sTmp, pSeqLim, pSeqExpo, lScanTimeBasic, lScanTimeSatsEtc, 0, SeqBuildBlockOptfs::SPIR_CALC_TYPE_EPI, pSBBOptfsPrep))
                    {
                        TRACE_PUT1(TC_ALWAYS, TF_SEQ,"%s: The calc SPIR time has failed, prob due to an error in the tickle pulse  ", ptModule);
                        return false;
                    }
                    
                    // Update SBB scan time including SPAIR
                    mySeqLoop.setlSBBScanTime(lScanTimeSatsEtc);

                    bSuccess = mySeqLoop.calcFillTimes ( &sTmp, pSeqLim, pSeqExpo );
                    
                    double dInc = pSeqLim->getTR()[0].getInc();
                    
                    // Avoid division by zero
                    if ( fabs(dInc) < 1.e-10 )
                    {
                        // Unexpected error
                        // => trace also if pSeqLim->isContextPrepForBinarySearch()
                        TRACE_PUT0(TC_INFO, TF_SEQ, "ERROR in calculateTRTIFillTimes(): dInc == 0. ");
                        return false;
                    }
                    
                    // ("..........NOTE!!!!! we assume something about TR-limit-handling in MrUILink !")
                    if (mySeqLoop.getlTRneeded() >   10000) 
                    {
                        dInc *= 10.0;
                    }
                    if (mySeqLoop.getlTRneeded() > 1000000) 
                    {
                        dInc *= 10.0;
                    }
                    
                    long lNeededTR_new = (long)(0.5+dInc*ceil(mySeqLoop.getlTRneeded()/dInc));
                    
                    *plNeededTR = max(lNeededTR_new, *plNeededTR);
                }
            }
            else
            {
                *plNeededTR = pMrProt->tr()[0];
            }
            
            if (mySeqLoop.getlTIneeded()!=0)
            {
                *plNeededTI = mySeqLoop.getlTIneeded();
            }
            else                             
            {
                *plNeededTI = pMrProt->ti()[0];
            }
            
            bSuccess = true;
        }
        else
        {
            // no chance
            //
            bSuccess = false;
        }
  }
  else
  {
      *plNeededTR = pMrProt->tr()[0];
      *plNeededTI = pMrProt->ti()[0];
  }
  
  // benpos SMS-CAIPI -->
	//put back to full number of slices 
	mySeqLoop.setSlicesToMeasure( pMrProt->sliceSeries().size()  );
    // <-- benpos SMS-CAIPI

  return bSuccess;
}


#ifndef VXWORKS
// ------------------------------------------------------------------------------
// Functions   : fSEQConvProt
// ------------------------------------------------------------------------------
//               
// Description : try to convert protocols from previous software versions
//
//
// Return      : SEQU_NORMAL for success
//               SEQU_ERROR  for error
//
// ------------------------------------------------------------------------------
#ifndef EP2D_DIFF
NLS_STATUS fSEQConvProt (const MrProt &rMrProtSrc, MrProt &rMrProtDst )
#else
NLS_STATUS fSEQConvProt (const MrProt & /*rMrProtSrc*/, MrProt & /*rMrProtDst */)
#endif
{
    static const char * const ptModule             = {"fSEQConvProt"};
    NLS_STATUS                lRet                 = SEQU_NORMAL;

    // charm MR_00348199: actively set coilCombineMode for protocol conversion
    //      automatic protocol conversion to COILCOMBINE_ADAPTIVE_COMBINE is here disabled
    #ifndef EP2D_DIFF
      // convert to new parameters if protocol version is earlier than first VB15
      if(rMrProtSrc.getConvFromVersion() < 21510000 )
      {
    rMrProtDst.coilCombineMode(SEQ::COILCOMBINE_SUM_OF_SQUARES);
      }
    #endif

    return lRet;
}
#endif  // ndef VXWORKS


/*[ Function ****************************************************************\
*
* Name        : fSEQInit
*               
* Description : Defines the hard limits for the Seq/Change dialog.
*               
* Return      : An NLS status code.
*
\****************************************************************************/

/*] END: */

NLS_STATUS fSEQInit
(
    SeqLim *pSeqLim
)
{
  #undef  DEBUG_ORIGIN
  #define DEBUG_ORIGIN 0x00000100
  
  static const char *ptModule = {"fSEQInit"};
  
  NLS_STATUS lStatus = SEQU__NORMAL;
  
  mTRrun;

// ggg
pSeqLim->disableSAFEConsistencyCheck();

#ifdef EP2D_DIFF
  #ifndef VXWORKS
    #include <process.h>
    if (theConfig.getIntDefault (ptModule, "WaitForDebugger", 0)) {
      TRACE_PUT2 (TC_INFO, TF_SEQ, "%s HALT: Attach debugger to PID %d and press <ret> to continue", 
		  ptModule, _getpid()) ;
      char foo[256] ;
      gets (foo) ;    
    }
  #endif
#endif

  //-------------------------------------------------------------------------------------
  // Init the pace feedback class
  //-------------------------------------------------------------------------------------
  #ifdef PACE3D
      // JR inserted PACE3d parameter (starting feedback in ICE)
      //pSeqLim->setBold3dPace(1, 0);	  
      pSeqLim->setBold3dPace(SEQ::ON, SEQ::OFF);	

      // init 3D PACE class
      gPaceFeedback.Init();
      
      #ifndef PASL   
      // handle time required for wakeup-eventblock within TR
      mySeqLoop.setAdditionalTimeForExtraEBinTRUsec( TWAKEUP );
      #else
      mySeqLoop.setAdditionalTimeForExtraEBinTRUsec( TWAKEUP + Pasl.getDurationPerRequest());
      #endif      
  #else
      pSeqLim->setBold3dPace(SEQ::OFF);
      #ifndef PASL  
      // default no plugin within TR
      mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(0);
      #else
      mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(Pasl.getDurationPerRequest());
      #endif
  #endif

  //-------------------------------------------------------------------------------------
  // sequence should not be executed on Magnetom OPEN
  //-------------------------------------------------------------------------------------
  pSeqLim->setNotSupportedSystemTypes("007");

  //-------------------------------------------------------------------------------------
  // sequence hint text
  //-------------------------------------------------------------------------------------
  {
    char t[512];
    sprintf(t,"%s (compile date: %s)",ptVariant,(char*)__DATE__);
    pSeqLim->setSequenceHintText(t);
  }

  //-------------------------------------------------------------------------------------
  // set sequence type
  //-------------------------------------------------------------------------------------
  pSeqLim->setSequenceType( SEQ::SEQUENCE_TYPE_EPI );
  pSeqLim->getSequenceType().setDisplayMode( SEQ::DM_OFF );

  //-------------------------------------------------------------------------------------
  // initialize maximum allowed data rate for online regridding
  //-------------------------------------------------------------------------------------
  dMaxAllowedDataRate = accMeasPerm.get_sICL_flAbsoluteMaxDataRateForOnlineRegriddingOnActualMRIR();

  //-------------------------------------------------------------------------------------
  // general settings
  //-------------------------------------------------------------------------------------
  pSeqLim->setMyOrigFilename                     (                    (char *) __FILE__);

  //-------------------------------------------------------------------------------------
  // change default gradient performance of excitation SBB
  //-------------------------------------------------------------------------------------
  double dRiseTimeMin = theGPA.getGradMinRiseTimeAbsolute();
  double adRiseTimes [3] = {dRiseTimeMin, dRiseTimeMin, dRiseTimeMin};
  
#ifndef VXWORKS
#pragma message ("NOTE: Workaround - modification of SBBBinomialPulses desired.")
#endif
  double dGradMaxMagnExt = theGPA.getGradMaxAmpl( SEQ::GRAD_FAST );
  double dGradMaxMagnBin = dGradMaxMagnExt;
  if (SysProperties::isTrio())
  {
      dGradMaxMagnBin = theGPA.getGradMaxAmplNominal();
  }
  double adMagnitudesExt[3] = {dGradMaxMagnExt, dGradMaxMagnExt, dGradMaxMagnExt };
  double adMagnitudesBin[3] = {dGradMaxMagnBin, dGradMaxMagnBin, dGradMaxMagnBin };

  WE_121.setMinRiseTimes (adRiseTimes    ,SBBBinomialPulses_GRAD_PERF_BINOMIAL   );
  WE_121.setMinRiseTimes (adRiseTimes    ,SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF);
  WE_121.setMaxMagnitudes(adMagnitudesBin,SBBBinomialPulses_GRAD_PERF_BINOMIAL   );
  WE_121.setMaxMagnitudes(adMagnitudesExt,SBBBinomialPulses_GRAD_PERF_EXTERNAL_RF);

  //-------------------------------------------------------------------------------------
  // protocol independent configuration of water excitation
  //-------------------------------------------------------------------------------------
  WE_121.setBandwidthTimeProduct (5.2); 

  #ifdef EP2D_SE
    // invert GS polarity of WE pulse compared to refocusing pulse to avoid 3rd arm artefact (checked by SeqUT)
    WE_121.setRequiredGSPolarity (-1.0);    
  #endif

  WE_121.setUsePossibleBandwidthTimeProduct(true);    

  //-------------------------------------------------------------------------------------
  // Set gradient mode options with default "fast"
  // The "normal" mode option has been introduced to allow a reduced slew-rate 
  // to be set for under-voltage situations (CHARM 333764) - see fSEQPrep.
  //-------------------------------------------------------------------------------------
  pSeqLim->setGradients (SEQ::GRAD_FAST, SEQ::GRAD_FAST_GSWD_RISETIME, SEQ::GRAD_NORMAL, SEQ::GRAD_NORMAL_GSWD_RISETIME);
  
  //-------------------------------------------------------------------------------------
  // EPI kernel configuration
  //-------------------------------------------------------------------------------------
  
  // set polarity of first RO pulse in echo-train to positive
  EPIKernel.setStartImagingReadOutWithNegativeGradient( false);

  // specify that three phase correction echoes are acquired after each excitation
  EPIKernel.setInternalPhaseCorrection( true );

  // activate RO ramp-sampling
  EPIKernel.setUseRegriddingForRO( true );

  // set slew rate to 95% of "fast" value specified in MeasPerm section (CHARM 32446)	
  // note that this can be changed to 80% by switching to "normal" gradient mode (see fSEQPrep)
  // Also note that this call also sets the default gradient performance (SBBEPIReadout)
  EPIKernel.setMinRiseTimeScaleFactor( 1.05 );

  //-------------------------------------------------------------------------------------
  // we do not use calculation limits
  // => full available gradient performance is used all the time
  //-------------------------------------------------------------------------------------
  myCalcLimits.resetAllLimits();

  if(!EPIKernel.setPointerToCalculationLimits(&myCalcLimits))
  {
      TRACE_PUT1_NLS(TC_INFO, TF_SEQ, "%s: EPIKernel.setPointerToCalculationLimits failed.", ptModule, EPIKernel.getNLSStatus());
      return EPIKernel.getNLSStatus();
  }

  //-------------------------------------------------------------------------------------
  // preparation of osc-bit
  //-------------------------------------------------------------------------------------
  if(!EPIKernel.setUseSyncBits(true)) // using default arguments from SeqBuildBlockEPIKernel::setUseSyncBits
  {
      TRACE_PUT1_NLS(TC_INFO, TF_SEQ, "%s: EPIKernel.setUseOscBit failed.", ptModule, EPIKernel.getNLSStatus());
      return EPIKernel.getNLSStatus();
  }

  //-------------------------------------------------------------------------------------
  // set call-back function for excitation SBB configuration
  //-------------------------------------------------------------------------------------
  if(!EPIKernel.setConfigureExcitationSBBFunction(&fConfigureExcitationForEPIKernel))
  {
      TRACE_PUT1_NLS(TC_INFO, TF_SEQ, "%s: EPIKernel.setConfigureExcitationSBBFunction failed.", ptModule, EPIKernel.getNLSStatus());
      return EPIKernel.getNLSStatus();
  }
  
  //-------------------------------------------------------------------------------------
  // set call-back functions for additional RTEB called by the kernel
  //-------------------------------------------------------------------------------------
  #ifdef EP2D_DIFF
    EPIKernel.setRTEBPlugIn (&fPrepDIFFPlugInForEPIKernel,&fRunDIFFPlugInForEPIKernel);
    if ( theConfig.getIntDefault(ptModule, "PrephaseAfterSBBDiffusion", 1) ) {
      // TRACE_PUT0 (TC_INFO, TF_SEQ, "EPIKernel.setPrephaseAfterRTEBPlugIn (true)" );      
      EPIKernel.setPrephaseROAfterRTEBPlugIn (true);
      EPIKernel.setPrephaseBlipsAfterRTEBPlugIn (true);
    }
    else 
      TRACE_PUT0 (TC_INFO, TF_SEQ, "EPIKernel.setPrephaseAfterRTEBPlugIn (false)" );      
    
    // Display master INI file  
    char tINIFileName[128], bas[128] ;
    theConfig.fGetMasterINIFilename(tINIFileName, bas);
    if (strlen(tINIFileName))
      TRACE_PUT2 (TC_INFO, TF_SEQ, "%s: iniFile=%s", ptModule, tINIFileName) ;
  #endif

  #ifdef EP2D_SE
    EPIKernel.setRTEBPlugIn(&fPrepSEPlugInForEPIKernel,&fRunSEPlugInForEPIKernel);
  #endif

  #ifdef PASL 
    EPIKernel.setRTEBPlugIn(&fPrepPaslPlugInForEPIKernel,&fRunPaslPlugInForEPIKernel);
    EPIKernel.setPrephaseBlipsAfterRTEBPlugIn (false); // benpos ASL true --> false
  #endif

  //-------------------------------------------------------------------------------------
  // protocol independent configuration of SeqLoop
  //-------------------------------------------------------------------------------------
  mySeqLoop.setPerformTRFill        (FALSE);
  #ifndef PASL   
  mySeqLoop.setPerformSATs          (TRUE);
  #else
  // we need to do the SATs ourselves
  mySeqLoop.setPerformSATs          (FALSE);
  #endif
  mySeqLoop.setEffectiveTRForR_CSat (FALSE,FALSE);
  mySeqLoop.setbSpoilGradAfterKernel(FALSE);
  mySeqLoop.setdDistFacMinIfConcYes (-1.0);
  mySeqLoop.setbHandleTRTIConflict  (TRUE); 
  mySeqLoop.initRegistryEntries     ();
  mySeqLoop.setTrigHaltSingleShot   (FALSE);
  mySeqLoop.setPerformOscBit        (FALSE);

  #if defined EP2D_SE || EP2D_DIFF
      mySeqLoop.setScaleIRThickness   (2.0);
  #endif

  //-------------------------------------------------------------------------------------
  // read debug-masks 
  //-------------------------------------------------------------------------------------
  #ifndef VXWORKS
    lDebug_SEQ_fSEQPrep = getMaskFromRegistry("DEBUG_SEQ_fSEQPrep");
  #endif

  //-------------------------------------------------------------------------------------
  // standard epi hard limits and solve handlers
  //-------------------------------------------------------------------------------------
  if (!fEPIStdInit(pSeqLim, (SeqBuildBlockEPIReadOut*)&EPIKernel, (ReorderInfo*)&REOInfo, &configInfo) )
  {
      textTR("fEPIStdInit failed");
      return SEQU_ERROR;
  }

  //-------------------------------------------------------------------------------------
  // standard single-shot epi hard limits
  //-------------------------------------------------------------------------------------
  pSeqLim->setEPIFactor                (       1,         512,SEQ::INC_SINGLESHOT,  128);
  pSeqLim->setMultiSliceMode           (                           SEQ::MSM_INTERLEAVED);
  pSeqLim->setDelayTimeInTR            (       0,    30000000,        1000,           0);
  pSeqLim->setSliceSeriesMode          (SEQ::INTERLEAVED,SEQ::ASCENDING,SEQ::DESCENDING);
  pSeqLim->setMultipleSeriesMode       (SEQ::MULTIPLE_SERIES_OFF, SEQ::MULTIPLE_SERIES_EACH_MEASUREMENT);
  pSeqLim->setSliceThickness           (   0.25,      10.000,       0.10,       5.000);  //benpos min 1->0.25
  
  // base resolution with flexible adjustment of incrememnt
  pSeqLim->setBaseResolution           (64, 512, SEQ::INC_NORMAL, 64                  );
  

  #ifdef SUPPORT_iPAT_a_ep2d
  {
    // set default SeqLim for PAT: PATMode, RefScanMode, AccelFactorPE, RefLinesPE
    fPATSetDefaultSeqLim(pSeqLim);
    // EPI uses the ExtraRefScanMode:
    pSeqLim->setRefScanMode    (SEQ::PAT_REF_SCAN_EXTRA);
  }
  #endif

  //-------------------------------------------------------------------------------------
  // variant specific limits
  //-------------------------------------------------------------------------------------

  // benpos WIP 
  // add special card parameters needed for multi-echo EPI

  // --> the starting position: BEGIN_PARAMETER_MAP(pSeqLim, first_alFreeParam, first_adFreeParam)
  BEGIN_PARAMETER_MAP(pSeqLim, 0, 3 );

#ifndef PASL  
	PARAM("save individual echoes ", &b_Checkbox_SaveSeparateEchoes, false, "Click here to also save the separate echoes.");
	PARAM_SELECT("echo combination", &l_SelectionBox_EchoCombineMode, 0, "Echoes can be combined online in different ways \nif only Magnitude reconstruction is selected:  \n\nBOLD weighted: \n  BOLD optimised weights using voxelwise weights given by the \n  product of echo time and sigal at that TE:\n  weight_n = S_n * TE_n (the weights are normalised to one) \n  optimal weights are obtained dynamically from each volume \n\nfixed weights: \n  constant weights used for all voxels \n  optimised for the given T2*: \n  weight_n=TE_n*exp(-TE_n/T2*) (normalised)  \n\nsimple average: \n  echoes are NOT weighted and simply averaged (use at 7T)  \n\nlinear with TE: \n  weighting proportional to echo time \n\nPlease refer to:\nPoser et al, MagnResonMed 2006 Jun;55(6):1227-35 \nPoser and Norris, Neuroimage 2009 May;45(4):1162-72 "); 
		OPTION("none"      ,0);	
		OPTION("BOLD weighted"        ,1);
		OPTION("fixed weights"        ,2);
		OPTION("simple average"       ,3);
		OPTION("linear with TE"       ,4);
	PARAM_SELECT_END ();
	PARAM("optimise weights for T2*", "ms", &l_LongEntry_T2starweight, 10, 100, 1, 30,  "gray matter T2* value for which to optimize \nat 3T use ~30-35ms \nat 1.5T use ~40-50ms \ndepending on your brain region of main interest") ;	
#else	
	PARAM_GROUP();
		PARAM("Ampl", "%", &l_LongEntry_trfociAmpl, 0, 501, 1, 100,  "Renzos tr foci Pulse parameter \nSpielplatz 4");	
		PARAM("BWDTH","% 3.1kHz", &l_LongEntry_trfociBWDTH, 10, 1000, 1, 300,  "Renzos tr foci Pulse parameter \nSpielplatz 5 \n relative to 3.1 kHz \nIch empfele 300 fuer ASL");
	//	PARAM("Offset", "mm", &d_trfociOFFSET, -100, 100, 0.1, 0,  "Renzos tr foci Pulse parameter \nSpielplatz 6 \ndo not change \nread only \nonly for debugging ");
		PARAM("thickness", "%", &l_LongEntry_trfociD, 20, 1004, 1, 100,  "Renzos tr foci Pulse parameter \nSpielplatz 7 \nRelative to 7cm");
    //PARAM_GROUP_END();


	//PARAM_GROUP();
	PARAM("Volumes per TI", " ", &l_LongEntry_VolumesPerInversion, 1, 25, 1,1,  " Multi-TI readout: number of imaging volummes per inversion");
	PARAM("FatSat flip angle", "deg", &l_LongEntry_FatSatFlipAngle, 0, 180, 1, 110,  "FatSat flip angle, default=110deg") ;	
	PARAM("RF pulse duration", "us", &l_LongEntry_RfDuration, 512, 10240, 512, 2560,  "RF pulse duration: Shorter pulses allow thinner slices \n(default 2560)");
	PARAM_GROUP_END();
	
	//ggg SKIP_PARAM();
	
#endif	
	
	
	PARAM_GROUP();
		PARAM("delay after echo 1", "ms", &l_LongEntry_delayBetweenFirstAndSecondEcho, 0, 100, 1,0,  "delay between first and sencond echo");
		PARAM("delay after echo 2...n", "ms", &l_LongEntry_delayBetweenRemainingEchoes, 0, 100, 1,0,  "delay between the other echoes");
	PARAM_GROUP_END();
	
	
	
	// benpos SMS-CAIPI -->
	PARAM_GROUP();
		PARAM("SMS factor", " ", &l_LongEntry_SmsFactor, 1, 16, 1,1,  "Slice acceleration factor");
		PARAM("CAIPI shift", " ", &l_LongEntry_SmsShiftFactor, 1, 8, 1,1,  "FoV shift between adjacent simultaneous slices is (1/shift)");
	PARAM_GROUP_END();
	
	PARAM_GROUP();
		PARAM("SMS online recon ", &b_Checkbox_SmsOnlineRecon, true, "Click here to enable online SMS reconstruction" ); 
		PARAM("SMS-RF phase optim.", &b_Checkbox_SmsPhaseOptimise, true,  "RF phase optimisation for the SMS pulses (see e.g. Wong ISMRM 2012)");
	PARAM_GROUP_END();
	// <-- benpos SMS-CAIPI

	
	PARAM("log physio files ", &b_Checkbox_LogPhysio, false, "Click here to log Pulse, Breathing and ExtTrig\nconnect the appropriate sensors \nlogfiles are located in C:/Medcom/log/Physlog \nand carry the timestamp of the start of the scan.");
	PARAM_GROUP();
		PARAM("altern z-shim", "uT/m", &l_LongEntry_AlternZshim, -450, 450, 1, 0,  "alternating z-shim gradient, alternates between volumes");
		PARAM("fixed z-shim", "uT/m", &l_LongEntry_FixedZshim, -450, 450, 1, 0,  "fixed z-shim gradient. ");
	PARAM_GROUP_END();
	
	PARAM_SELECT("EPI phase correction", &l_SelectionBox_PhaseCorrMode, 0, " type of EPI phase correction \n'local' tends to achieve better ghost correction \nespecially when using the 32ch coil and in case of frequency drifts \n 'pixel-by-pixel' does simple point-by-point phase conjugation ' ") ;
		OPTION("normal"     ,0);	
		OPTION("local"      ,1);
		OPTION("pxl-by-pxl" ,2);
	PARAM_SELECT_END ();

	PARAM_SELECT("PAT refscan mode", &l_SelectionBox_RefScanMode, 0, "GRAPPA reference lines are acquired with an EPI readout in either a segmented or a single-shot fashion \nsince the same EPI readout is used for reference lines as for imaging scans \nthis means tha much fewer reference lines can be acquired in single-shot mode \nbut it tends to be less vulnerable to ghosting, especially  if there is motion during the reference scans... \n\nuse of a Flash refence scan may leave more residual ghost \nbut the ghost tends to be temporally more stable  ") ;
		OPTION("segmented"    ,0);	
		OPTION("single-shot"  ,1);
		OPTION("FLEET"  ,2);
		OPTION("Flash"  ,3);
	PARAM_SELECT_END ();


	PARAM_GROUP();
		PARAM("FLEET dummies", " ", &l_LongEntry_FleetDummies, 1, 15,  1 , 5,  "FLEET dummy scans");
		PARAM("FLEET flip angle", " ", &l_LongEntry_FleetFlipAngle, 1, 90,  1 , 15,  "flip angle used for FLEET reference scans");
//must be long		PARAM("Flash refscan BaseRes", " ", &d_DoubleEntry_SSBPatRefscanBaseRes, 32, 256,  2 , 64,  "sets BaseResolution for Flash based reference scan");
//long! 		PARAM("Flash refscan TE", "ms ", &d_DoubleEntry_SSBPatRefscanTE, 3, 20, 1, 4,  "sets TE for Flash based reference scan");
//long!		PARAM("Flash refscan BW ", "Hz/px ", &d_DoubleEntry_SSBPatRefscanBW, 100, 2000, 10,300,  "sets bandwidth for Flash based reference scan");
	PARAM_GROUP_END();

	
	//PARAM_GROUP();
		// ggg moved PARAM("RF pulse duration", "us", &l_LongEntry_RfDuration, 512, 10240, 512, 2560,  "RF pulse duration: Shorter pulses allow thinner slices \n(default 2560)");
	//	PARAM("RF BWTP", "", &d_DoubleEntry_RfBWTP, 1, 16, 1, 5,  "RF bandwidth-time product (default 5.2)");
	//PARAM_GROUP_END();
	
	PARAM("FFT scale", " ", &d_DoubleEntry_FFTscale, 0, 10, 0.1,1,  "additional FFT scale factor");

	// ggg
	
	//PARAM_GROUP();
	PARAM("RO samples", 	"",  &paramLongROSamples, 1024, 	81920, 1024,	7168,	"");
	//PARAM("Interleaves",	"",  &paramLongInterleaves, 1,	 122,	1,			1,		"");
	//PARAM("Spiral Slew Rate","T/m/s", &paramLongSpSlewRate, 	1,	  MAX_SLEW_RATE_SP,	1, 125, "");
	//PARAM_GROUP_END();
	//VERIFY_OFF(&paramLongROSamples);
	//VERIFY_OFF(&paramLongInterleaves);
	/*
	PARAM("Spiral Peak Grad","mT/m",	 &paramLongSpGradAmp,	1,	 MAX_GRAD_AMP_SP, 1,   25,	"");
	PARAM("Spiral Type",		"", 	&paramLongSpiralType,	0,		 3, 	1,		0,	 "");
	PARAM("Spiral BW",		"", 	&paramLongSpiralBW,	200000,		 1000000, 	1,		400000,	 "");

	PARAM("Flat First ZC",		"", 	&paramLongZCToFlat,	0,		 10, 	1,		0,	 "");*/
	
	
	
	//PARAM_GROUP_END();
	/*
	VERIFY_OFF(&paramLongROSamples);
	VERIFY_OFF(&paramLongInterleaves);
	VERIFY_OFF(&paramLongSpSlewRate);
	VERIFY_OFF(&paramLongSpGradAmp);
	VERIFY_OFF(&paramLongSpiralType);
	VERIFY_OFF(&paramLongSpiralBW);
	VERIFY_OFF(&paramLongZCToFlat);
	*/
	PARAM_GROUP();

	PARAM("VD",	"", 	&paramDoubleVD	,	-1.0,	 10.0,		.1,		-1.0,	"");
	PARAM("Undersampling fac",	"", 	&paramDoubleAccR	,	0,	 20.0,		.1,		1.0,	"");
	PARAM("CAIPI shift (mm)",	"", 	&paramDoubleEchoShift	,	-50,	 50.0,		.1,		0.0,	"");
	PARAM("CAIPI period (us)",	"", 	&paramDoubleCAIPIPeriod	,	0,	 9000.0,		.1,   300.0,	"");
	PARAM("CAIPI delay (us)",	"", 	&paramDoubleCAIPIDelay	,	0,	 9000.0,		.01,		0.0,	"");
	
	
	PARAM_GROUP_END();
	
	VERIFY_OFF(&paramDoubleVD);
	VERIFY_OFF(&paramDoubleAccR);
	VERIFY_OFF(&paramDoubleEchoShift);
	VERIFY_OFF(&paramDoubleCAIPIPeriod);
	

	PARAM_GROUP();
	PARAM("Interleaves",	"",  &paramLongInterleaves, 1,	 122,	1,			1,		"");
	PARAM("Spiral Slew Rate","T/m/s", &paramLongSpSlewRate, 	1,	  MAX_SLEW_RATE_SP,	1, 125, "");
	PARAM("Spiral Peak Grad","mT/m",	 &paramLongSpGradAmp,	1,	 MAX_GRAD_AMP_SP, 1,   25,	"");
	PARAM("Spiral Type",		"", 	&paramLongSpiralType,	-3,		 24, 	1,		0,	 "");
	PARAM("Spiral BW",		"", 	&paramLongSpiralBW,	0	,		 3000000, 	1,		400000,	 "");
	PARAM("Flat First ZC",		"", 	&paramLongZCToFlat,	0,		 50, 	1,		0,	 "");
	PARAM_GROUP_END();
	VERIFY_OFF(&paramLongROSamples);
	VERIFY_OFF(&paramLongSpSlewRate);
	VERIFY_OFF(&paramLongSpGradAmp);
	VERIFY_OFF(&paramLongSpiralType);
	VERIFY_OFF(&paramLongSpiralBW);
	VERIFY_OFF(&paramLongZCToFlat);

	TOOLTIP_CUSTOM(&paramLongROSamples, _WIP_LONG_GetToolTipVD);
	//TOOLTIP_CUSTOM(&paramDoubleVD, _WIP_DOUBLE_GetToolTipVD);*

	//PARAM_GROUP_END();
	//VERIFY_OFF(&paramLongROSamples);
	//VERIFY_OFF(&paramLongInterleaves);
	/*
	PARAM("Spiral Peak Grad","mT/m",	 &paramLongSpGradAmp,	1,	 MAX_GRAD_AMP_SP, 1,   25,	"");
	PARAM("Spiral Type",		"", 	&paramLongSpiralType,	0,		 3, 	1,		0,	 "");
	PARAM("Spiral BW",		"", 	&paramLongSpiralBW,	200000,		 1000000, 	1,		400000,	 "");

	PARAM("Flat First ZC",		"", 	&paramLongZCToFlat,	0,		 10, 	1,		0,	 "");*/
	
	 
  END_PARAMETER_MAP;

  gSpiral_Init(pSeqLim);

  // ggg
  pSeqLim->disableSAFEConsistencyCheck();

  //make the Contrast variable editable, and set limit to 15 echoes
  //this overwrite some defaults set in epi_StdUILink.cpp
  int max_contrasts = 15;
  pSeqLim->getContrasts().setDisplayMode(SEQ::DM_EDIT);
  pSeqLim->setContrasts( 1, max_contrasts ,  1,  1) ;
  for (int i = 1; i < max_contrasts; i++)
  {
	   pSeqLim->setTE(i,1000, 2500000,100, 1000);
  }

// benpos do something like this to make the phaseFOV larger than the readFOV
//pSeqLim->setReadoutFOV                ( 20, 500, 0.1, 400 ); 
//pSeqLim->setPhaseFOV                  ( 20, 1000, 0.1, 600 );
//



  //-------------------------------------------------------------------------------------
  #ifdef EP2D_DIFF
  //-------------------------------------------------------------------------------------
   pSeqLim->setSequenceHintText( (char *) "\n\
      Sequence: ep2d_diff "
#ifdef WIP
"WIP sequence. Not for clinical use." 
#endif
"\n\
   Application: EPI-based diffusion weighted imaging and DTI.\n\
        Basics: 2D single shot EPI; asymmetric sampling in phase direction;\n\
                Diffusion encoding with four bipolar gradients and double spin echo;\n\
                pixelwise phase correction using three measured projection lines;\n\
                read-out module optimized for gradient performance\n\
         Build: "__DATE__"  "__TIME__"\n" );

    pSeqLim->setPhasePartialFourierFactor (SEQ::PF_6_8, SEQ::PF_7_8, SEQ::PF_OFF, SEQ::PF_5_8) ;
    pSeqLim->setInversion                 (SEQ::INVERSION_OFF,SEQ::SLICE_SELECTIVE);   /* FLAIR ;-) */
    pSeqLim->setTI                     (0,       0,     9000000,         100,     2500000);
    pSeqLim->setRFPulseType               (SEQ::RF_NORMAL) ;
    pSeqLim->setFlipAngle                 (90.0, 90.0, 1.0, 90.0) ;
    pSeqLim->setNoiseLevel                (0, 4095, 1, 20);
    if (theMSU.getNominalB0() > LOW_FIELD_B0)
    {
        pSeqLim->setFatSuppression       (SEQ::FAT_SATURATION, SEQ::WATER_EXCITATION, SEQ::FAT_SUPPRESSION_OPTIMAL, SEQ::FAT_SUPPRESSION_OFF);
    }

    // invalid protocol if distortion correction is active
    pSeqLim->setDistortionCorrMode     (SEQ::DISTCORR_NDIS,SEQ::DISTCORR_DIS2D); 




    // modify maximum TE and TR after introduction of modified PAT selection solve-handler (CHARM 315966)
    // previous values did not allow solve-handler to work correctly with large number of slices
    pSeqLim->setTE                     (0,    1000,      400000,         100,      200000);
    pSeqLim->setTR                     (0,   10000,    30000000,        1000,      500000);

    // Do not display the following UI parameters:
    pSeqLim->getFlipAngle().setDisplayMode (SEQ::DM_OFF);
    // pSeqLim->setMultipleSeriesMode(SEQ::MULTIPLE_SERIES_EACH_MEASUREMENT);
    // pSeqLim->getMultipleSeriesMode().setDisplayMode(SEQ::DM_OFF);


    // Initialize SBBDiffusion
    char tError[512];
    char tIniFileName[512];
    if ( theConfig.fGetMasterINIFilename(tIniFileName, tError) ) {
      TRACE_PUT2 (TC_INFO, TF_SEQ, "%s WARNING: This sequences uses an ini file '%s'", ptModule, tIniFileName) ;
      Diff.setIniAccess(&theConfig);
      // If the user has specified an ini file, SIEMENS warranty is void:
      pSeqLim->setSequenceOwner("USER");
    }    
    if (! Diff.init(pSeqLim) ) {
      // Diff does not set NLS-Status, but returns always true. So we should never get into here ;-)
      TRACE_PUT1(TC_INFO, TF_SEQ, "%s ERROR: Diff.init failed",ptModule);
      return SEQU_ERROR;          
    }
  
    //file containing the default postprocessing protocol (EVAProtocol)
    #ifndef VXWORKS
     pSeqLim->setDefaultEVAProt (_T("%SiemensEvaDefProt%\\DTI\\DTI.evp"));
    #endif

  //-------------------------------------------------------------------------------------
  #endif

  //-------------------------------------------------------------------------------------
  #ifdef EP2D_FID
  //-------------------------------------------------------------------------------------
    pSeqLim->setTD                     (0,       0,    10000000,         100,           0);
    pSeqLim->setPhasePartialFourierFactor(SEQ::PF_OFF,SEQ::PF_7_8,SEQ::PF_6_8, SEQ::PF_5_8, SEQ::PF_HALF     ); //benpos added 5/8 and 4/8 (half) for those who dare...
    pSeqLim->setInversion                (                             SEQ::INVERSION_OFF);
    pSeqLim->setRFPulseType              (                                 SEQ::RF_NORMAL);
    pSeqLim->setRFSpoiling                                                      (SEQ::OFF);
    pSeqLim->getRFSpoiling().setDisplayMode                                  (SEQ::DM_OFF);

	// default number of repetitions must be consistent with EVA protocol (CHARM 339267)
    pSeqLim->setRepetitions              (       0,        4095,           1,          14);  

    #ifdef BOLD
      pSeqLim->setPELines                (      32,         256,           1,          64); //benpos max 128 --> 256
      pSeqLim->set2DInterpolation        (SEQ::OFF                                       );
      pSeqLim->setRepetitions            (       0,        8191,           1,          0 ); //benpos max 4095 --> 8191 and default 19 --> 0

	  pSeqLim->setFatSuppression         ( SEQ::FAT_SATURATION, SEQ::FAT_SUPPRESSION_OFF); //SEQ::WATER_EXCITATION); //benpos caipi (cannot multi-plex WE121 pulses! )
	  
      // mosaic image format does not support multi-slice-multi-angle
      pSeqLim->disableMSMA               (                                               );

      // long-term averaging is not supported by IceProgramOnline2D (CHARM 308391)
      pSeqLim->setAverages               (       1,           1,           1,           1);

	  // set flag to exclude SET and REP dimensions from data amount calculation in Sequence.cpp this is appropriate
	  // for ep2d_bold and ep2d_pace sequences, which perform inline image calculation without raw data storage
	  pSeqLim->setInteractiveRealtime (SEQ::ON);

      // invalid protocol if distortion correction is active
      pSeqLim->setDistortionCorrMode     (SEQ::DISTCORR_NDIS); 

      //file containing the default postprocessing protocol (EVAProtocol)
      #ifndef VXWORKS
        //benpos WIP
		// the modified evp currently has to into the CustomerSequence directory - a save and regularly back-up place to put it
	  pSeqLim->setDefaultEVAProt (_T("%SiemensEvaDefProt%\\BOLD\\t-test_10B10A_moco.evp"));
 //xxxxxxxx		pSeqLim->setDefaultEVAProt (_T("%CustomerSeq%\\BP_evaprot_WIP.evp"));
		

      #endif

      // sequence ep2d_pace can only be run with fatsat on (CHARM 309540)
      #ifdef PACE3D
         pSeqLim->setFatSuppression      ( SEQ::FAT_SATURATION );
      #endif // PACE3D

    #endif // BOLD

    #ifdef PASL  

      pSeqLim->setInteractiveRealtime (SEQ::ON); //benpos 
      pSeqLim->setBaseResolution         (      64,         512, SEQ::INC_NORMAL,      64); // enabled for pasl
      pSeqLim->setPELines                (      32,         512,           1,          64);
      pSeqLim->set2DInterpolation        (SEQ::OFF                                       );
      pSeqLim->setRepetitions            (       0,        4*4095,           1,          20); // should be only odd !! no! even too!
      pSeqLim->setSliceSeriesMode        (SEQ::ASCENDING,SEQ::DESCENDING,SEQ::INTERLEAVED);
      pSeqLim->setTR                     (0,   10000,    30000000,         1,     5000000); // pushed up for pasl 
      pSeqLim->setIntro                  (SEQ::ON, SEQ::OFF                              ); // for pace
      pSeqLim->setDelayTimeInTR          (       0,           0,        1000,           0);  // disabled
      
      // mosaic image format does not support multi-slice-multi-angle
      pSeqLim->disableMSMA               (                                               );

      // long-term averaging is not supported by IceProgramOnline2D (CHARM 308391)
      pSeqLim->setAverages               (       1,           1,           1,           1);

      // sequence ep2d_pace can only be run with fatsat on (CHARM 309540)  
      pSeqLim->setFatSuppression         ( SEQ::FAT_SATURATION, SEQ::FAT_SUPPRESSION_OFF );    //SEQ::FAT_SUPPRESSION_OFF, SEQ::WATER_EXCITATION
      
	   // allow two modes of fat saturation, weak - only one fat sat before all slices, strong - fat sat every slice
      pSeqLim->setFatSatMode             ( SEQ::FAT_SAT_WEAK, SEQ::FAT_SAT_STRONG ); // (CHARM 366331)

      // set Rsat for display of labelling region
      pSeqLim->setRSats                  (       0,           2,           1,           0);
      pSeqLim->setRSatThickness          (  100.000,    200.000,       1.000,     100.000);
      pSeqLim->setPSatMode               ( SEQ::PSAT_SINGLE_REG,     SEQ::PSAT_DOUBLE_REG);
      pSeqLim->setPSatThickness          (  100.000,    200.000,       1.000,     100.000);
      pSeqLim->setPSatGapToSlice         (   0.000,      50.000,       0.100,      25.000);  
      pSeqLim->setMTC                    (                                       SEQ::OFF);

      // set here mode for default protocol and starting values
      pSeqLim->setInversion              (                             SEQ::INVERSION_OFF); 
      pSeqLim->setAslMode                ( SEQ::ASL_FAIRQUIPSSII                         );       
      pSeqLim->setAslT2                  (       0,     4000000,         100,     900000); // ti()[0]
      pSeqLim->setAslT1                  (       0,     2000000,       25000,      700000); // ti()[1]
      pSeqLim->setAslT1Stop              (       0,     2000000,       25000,     1600000); // ti()[2]    
      pSeqLim->setAslFlowLimit (0.1, CRUSHGRAD_MAX_VELOCITY_ENC, 0.1, CRUSHGRAD_MAX_VELOCITY_ENC);

      //-------------------------------------------------------------------------------------
      // unclutter the UI, switching unused parameters off
      //-------------------------------------------------------------------------------------
      pSeqLim->getMTC().setDisplayMode(SEQ::DM_OFF);

      //file containing the default postprocessing protocol (EVAProtocol)
		#ifndef VXWORKS
        // pSeqLim->setDefaultEVAProt (_T("%SiemensEvaDefProt%\\ASL\\PASL_moco.evp")); // dimo - changed in order to match the previous sequence and debug reco artefact issues    
// DIMO		pSeqLim->setDefaultEVAProt (_T("%CustomerSeq%\\BP_evaprot_WIP.evp")); 
		 pSeqLim->setDefaultEVAProt (_T("%SiemensEvaDefProt%\\BOLD\\t-test_10B10A_moco.evp"));
 		#endif      
    #endif // PASL

    #ifdef PERF

      //file containing the default postprocessing protocol (EVAProtocol)
      #ifndef VXWORKS
         pSeqLim->setDefaultEVAProt (_T("%SiemensEvaDefProt%\\Perfusion\\GBP_PBP_TTP_relCBV.evp"));
      #endif

    #endif // PERF

  //-------------------------------------------------------------------------------------
  #endif // EP2d_FID

  //-------------------------------------------------------------------------------------
  #ifdef EP2D_SE
  //-------------------------------------------------------------------------------------
    pSeqLim->setTD                     (0,       0,    10000000,         100,           0);
    pSeqLim->setInversion                (        SEQ::INVERSION_OFF,SEQ::SLICE_SELECTIVE);
    pSeqLim->setTI                     (0,   10000,     5000000,         100,     2500000);
    pSeqLim->setRFPulseType              (                                 SEQ::RF_NORMAL);
    pSeqLim->setFlipAngle                (  90.000,      90.000,       5.000,      90.000);
    pSeqLim->getFlipAngle().setDisplayMode(SEQ::DM_OFF);

    // remove option for partial Fourier 4/8 due to poor image quality (CHARM 310630)
    pSeqLim->setPhasePartialFourierFactor(SEQ::PF_OFF,SEQ::PF_5_8,SEQ::PF_6_8,SEQ::PF_7_8);

    // restrict number of measurements for ep2d_se due to memory limitation in IceChest (CHARM 309396)
    pSeqLim->setRepetitions              (       0,         127,           1,           0);  

    //file containing the default postprocessing protocol (EVAProtocol)
    #ifndef VXWORKS
      pSeqLim->setDefaultEVAProt (_T("%SiemensEvaDefProt%\\Inline\\Inline.evp"));
    #endif
  //-------------------------------------------------------------------------------------
  #endif

  
	  
  //-------------------------------------------------------------------------------------	  
  // benpos WIP 
  // some additional changes:
  // - allow magnitude AND phase data
  // - set sequence ownder to Donders Institute
  pSeqLim->setReconstructionMode (SEQ::RECONMODE_MAGNITUDE,  SEQ::RECONMODE_MAGN_PHASE );
  pSeqLim->setSequenceOwner("Poser@MBIC, Maastricht, NL");
  //-------------------------------------------------------------------------------------
	  
	  
  //-------------------------------------------------------------------------------------
  // register solve-handeler IR-Pulse, TI getLimit- and solve-handler
  //-------------------------------------------------------------------------------------
  #if (defined EP2D_DIFF) || (defined EP2D_SE)
    if (!fEPIStdRegisterTIHandlers(pSeqLim, &calculateTRTIFillTimes))
    {
        textTR("fEPIStdRegisterTIHandlers failed");
        return SEQU_ERROR;
    }
  #endif

  //-------------------------------------------------------------------------------------
  // single shot epi sequences have problems with rephasing first PE-moment
  //-------------------------------------------------------------------------------------
  #ifndef VXWORKS
    SeqUT.DisableTestCase(lGpFirstMomentNotRephasedErr, RTEB_ORIGIN_fSEQRunKernel,
                          "single shot sequence, k-space center line not in first measured segment" );

	//benpos WIP - for z-shimmming also need to remove two checks on slice axis 
	SeqUT.DisableTestCase(lGsMomentNotRephasedErr, RTEB_ORIGIN_fSEQRunKernel,
                          "single shot sequence, k-space center line not in first measured segment" );
	SeqUT.DisableTestCase(lAddDimKzMomentCompareErr, RTEB_ORIGIN_fSEQRunKernel,
                          "single shot sequence, k-space center line not in first measured segment" );

  #endif

  
#ifdef SUPPORT_PACE
  //-------------------------------------------------------------------------------------
  // Initialization of SLFB
  //-------------------------------------------------------------------------------------
    if( !mySeqLoop.init(pSeqLim) )
    {
        TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d).",__FILE__,__LINE__);
        return mySeqLoop.getNLSStatus();
    }
#endif
  //-------------------------------------------------------------------------------------
  // register advanced handlers for the UI
  //-------------------------------------------------------------------------------------
#ifndef VXWORKS

  #ifdef EP2D_DIFF
    if( !fEPIRegisterEP2D_DIFFHandlers(pSeqLim) )
    {
        TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d).",__FILE__,__LINE__);
        return SEQU_ERROR;        
    }
  #endif

    //-----------------------------------------------------------------------------------
    // do not display 'Pause after Measurement'
    //-----------------------------------------------------------------------------------
    if(MrUILinkArray* _pPause = _search< MrUILinkArray > (pSeqLim, MR_TAG_MEASUREMENT_DELAY_TIMES))
    {
        _pPause->unregister();
    }

    //-----------------------------------------------------------------------------------
    // register solve handler for measurements:
    //
    // This solve handler became necessary due to CHARM 305893 (TR should always be 
    // physically correct). The sequence has to take into account the time for SUBFINI
    // and SUBSTRT between two repetitions for TR and measuerement time calculations.
    // If this is done only for multiple measurements, then TR may have to be increased,
    // if the number of measurements is increased from 1 to any other value.
    //
    // Unfortuneatly the UI already provides a solve-handler for measurements dealing with
    // problems of the EVA-protocols. Incorporating an original solve-handler into an
    // overloaded one is no trivial task, because reformatting of the message box is
    // necessary.
    //
    // So the problem is solved within the method SeqLoopEP2D::calcFillTimesOnly  by taking
    // into account the 340us for SUBFINI and SUBSTRT also for one measurement only.
    // Another possibility (which would be more transparent for the user) for solving the
    // problem would be to set the minimum for delayTimeInTR to 500us.
    //
    // SO THE SOLVE HANDLER IS NOT REGISTERED!
    //
    //-----------------------------------------------------------------------------------
    //if(LINK_LONG_TYPE *_measurements = _search<LINK_LONG_TYPE> (pSeqLim, MR_TAG_MEASUREMENTS))
    //{
    //    _measurements->registerSolveHandler (fEPIStdSolveLongParamConflict); 
    //}

    //-----------------------------------------------------------------------------------
    // register solve handler for 'delay in TR'
    //-----------------------------------------------------------------------------------
    if(LINK_DOUBLE_TYPE *_delayInTR = _search<LINK_DOUBLE_TYPE> (pSeqLim, MR_TAG_DELAY_IN_TR))
    {
        _delayInTR->registerSolveHandler (fEPIStdSolveDoubleParamConflict);
    }

  #ifdef PASL       
    //-------------------------------------------------------------------------------------
    // register the PASL UI parameters
    //-------------------------------------------------------------------------------------
    _PaslUI_Register(pSeqLim);

	//enable solve handling for fat sat mode and TR conflict (CHARM 366331)
	MrUILinkSelection<unsigned> *_fatsatMode       = _search< MrUILinkSelection<unsigned> >(pSeqLim, MR_TAG_FAT_SAT_MODE         );

	if (_fatsatMode) _fatsatMode->registerSolveHandler(fEPIStdSolveSelectionConflict);

  #endif

  #ifdef SUPPORT_iPAT_a_ep2d
    //-----------------------------------------------------------------------------------
    // register standard solve handlers for PAT
    //
    // CHARM 315966: restored PAT-specific selection solve-handler, which was not
    //               correctly copied from VA21B to VA25A archive.
    //-----------------------------------------------------------------------------------
    MrUILinkSelection<unsigned> *_PATMode                = _search< MrUILinkSelection<unsigned> >(pSeqLim, MR_TAG_PAT_MODE    );
    MrUILinkLimited<long>       *_PATAccelerartionFactor = _search< MrUILinkLimited<long>       >(pSeqLim, MR_TAG_PAT_ACC_PE  );
    MrUILinkLimited<long>       *_PATReferenceLines      = _search< MrUILinkLimited<long>       >(pSeqLim, MR_TAG_PAT_LINES_PE);

    if (_PATMode               ) _PATMode               ->registerSolveHandler(fEPIPATModeSolveSelectionConflict);
    if (_PATAccelerartionFactor) _PATAccelerartionFactor->registerSolveHandler(fEPIStdSolveLongParamConflict);
    if (_PATReferenceLines     ) _PATReferenceLines     ->registerSolveHandler(fEPIStdSolveLongParamConflict);
  #endif

	//-----------------------------------------------------------------------------------
    // register additional optional UI handlers (see epi_StdUILink.cpp)
    //-----------------------------------------------------------------------------------
	fEPIRegisterOptionalHandlers( pSeqLim );

#endif // of #ifndef VXWORKS
       
  //-------------------------------------------------------------------------------------
  // finished
  //-------------------------------------------------------------------------------------
  mTRend;
  return lStatus;
}


/*[ Function ****************************************************************\
*
* Name        : fSEQPrep
*               
* Description : Prepares everything that the sequence needs at run time.
*               
* Return      : An NLS status code.
*
\****************************************************************************/

/*] END: */

NLS_STATUS fSEQPrep
(
  MrProt   *pMrProt, // IMP: Measurement protocol
  SeqLim   *pSeqLim, // IMP: Sequence limits
  SeqExpo  *pSeqExpo // EXP: Returned values
)
{
  #undef  DEBUG_ORIGIN
  #define DEBUG_ORIGIN 0x00000200
  
  static const char *ptModule = {"fSEQPrep"};
  mTRrun;

  NLS_STATUS   lStatus                       = SEQU__NORMAL; // default return status
  double       dEnergy                       = 0.0;
  double       dBandwidthPerPixelPE          = 0.0;
  long         lEffectiveEchoSpacing         = 0;
  long         lDefaultNoOfRXChannels        = 256;
  long         lNeededTE                     =  -1;
  long         lNeededTI                     =  -1;
  long         lNeededTR                     =  -1;
  long         lRequiredPrepScans            = 0;
  long         lRepetition                   = 0;
  bool         bSuccess                      = 0;
  
  
  // benpos SMS-CAIPI -->
  long         lFullTRMeas                   = 0;
  long         lTotalMeas                    = 0;
  long         lReducedTRMeas                = 0;
  long         lReducedTRMeas_SecondMeas     = 0;
  double       dEnergySeqLoopPerSlice        = 0.0;
  long         lTR_SliceAccPrepMeas          = 0;
  //long         lCoolPauseExplicit            = 0;
  long         lScanTimeSatsEtc              = 0;
  long         lScanTimeBasic                = 0;
  // <-- benpos SMS-CAIPI
  
  //benpos WIP
  //get the sequence special card into the protocol
  PREPARE_PARAMETER_MAP(pMrProt, pSeqLim);

	//ggg
  gSpiral_ResetTRElements();

  #ifdef PASL    
  // ---------------------------------------------------------------------------
  // initialise UI parameters on the sequence/special card for PASL
  // and force error if WIP UI parameters are not yet initialised
  // ---------------------------------------------------------------------------
  #ifndef VXWORKS
    if( !_PaslUI_Init (pMrProt, pSeqLim) ) 
      return SEQU_ERROR;
  #endif
  #endif

  // ---------------------------------------------------------------------------
  // check that dwell time is not less than hardware limit of 1000ns
  // ---------------------------------------------------------------------------
  for ( long lI = 0 ; lI < pMrProt->contrasts() ; lI++ )
  {
    if ( pMrProt->rxSpec().realDwellTime()[lI] < 1000 )
    {
        if (!pSeqLim->isContextPrepForBinarySearch()) textTR("hardware sampling-rate limit exceeded");
        lStatus=SEQU_ERROR;
        goto FINISHED;
    }
  }

  // Prohibit simultaneous use of IR and SPAIR
  if ((pMrProt->preparationPulses().fatSuppression()==SEQ::FAT_SUPPRESSION_OPTIMAL) && !(pMrProt->preparationPulses().inversion()==SEQ::INVERSION_OFF)) 
  {
      lStatus = SEQU_ERROR;
      goto FINISHED;
  }
  
  

  //-------------------------------------------------------------------------------------
  // Set image format in epiConfig class (to either 'standard' or 'mosaic')
  //
  // For ep2d_bold and ep2d_pace sequences image format is always 'mosaic'.
  // 
  // For ep2d_fid and ep2d_se sequences image format is always 'standard'
  //
  // For ep2d_diff 'mosaic' image format is used for DTI scans 
  // (i.e. in 'MDDW' and 'Free' modes) and 'standard' format 
  // is used in the other modes.
  //
  // This information is accessed by UI handlers to control availability of distortion 
  // correction, which is not compatible with mosaic image format in ICE (VB13A).
  //
  // NOTE: this setting does not actually control the image format used during 
  //       reconstruction - it is just provides information for UI functions 
  //       in epi_StdUILink.cpp
  //-------------------------------------------------------------------------------------
  #ifdef BOLD

	  configInfo.setImageFormat( IMAGE_FORMAT_MOSAIC );

  #else

      if (   ( pMrProt->diffusion().mode() == SEQ::DIFFMODE_TENSOR )
		  || ( pMrProt->diffusion().mode() == SEQ::DIFFMODE_FREE   )
         )
	  {	  
		  // invalid protocol if distortion correction is active

		  if ( pMrProt->distortionCorrFilter().mode() != SEQ::DISTCORR_NDIS )
		  {
			  lStatus = SEQU_ERROR;
	
			  if (! pSeqLim->isContextPrepForBinarySearch() )
			  {
				  TRACE_PUT1(TC_INFO, TF_SEQ, "%s: Distortion correction not supported with MDDW or free diffusion mode ", ptModule );
			  }

			  goto FINISHED;
		  }

		  // set image format to 'mosaic'

		  configInfo.setImageFormat( IMAGE_FORMAT_MOSAIC );
	  }

	  else
	  {
		  configInfo.setImageFormat( IMAGE_FORMAT_STANDARD );
	  }

  #endif

  // ---------------------------------------------------------------------------
  // reset solve handler control
  // ---------------------------------------------------------------------------
  fEPIStdResetSolveHandlerControlTETITR();

  // ---------------------------------------------------------------------------
  // set default value for maximum numnber of receiver channels
  // ---------------------------------------------------------------------------
  pSeqExpo->setMaxReceiverChannels(lDefaultNoOfRXChannels);

  // ---------------------------------------------------------------------------
  // Calculate the rotation matrices and offsets for slices
  // ---------------------------------------------------------------------------
  lStatus = fSUPrepSlicePosArray (pMrProt, pSeqLim, asSLC);

  // error can only be caused by programming error
  // => trace also if pSeqLim->isContextPrepForBinarySearch()
  CheckStatusPB(lStatus,"fSUPrepSlicePosArray");

  // ---------------------------------------------------------------------------
  // Specify type of reference scan for PAT.
  //
  // For PAT factors greater than 2 a segmented reference scan is used. The 
  // primary reason for this is that the same timing is currently used for the 
  // PAT reference scans as for the imaging scans. Consequently with a high PAT 
  // factor there is a restriction on the number of reference lines that can be 
  // acquired with a single-shot reference scan. 
  // ---------------------------------------------------------------------------
  #ifdef SUPPORT_iPAT_a_ep2d
  {
	  if (    ( pMrProt->PAT().PATMode() != SEQ::PAT_MODE_NONE )
           && ( pMrProt->PAT().AccelFactPE() > 2               )
		 )
	  {	  
		  bSegmentedRefLines = true;
	  }
		  
	  else  
	  {	  
		  bSegmentedRefLines = false;
	  }
  }
  
  //benpos - override this normal behaviour 
  if (l_SelectionBox_RefScanMode == 0 )  bSegmentedRefLines = true;
  if (l_SelectionBox_RefScanMode == 1 )  bSegmentedRefLines = false;
  if (l_SelectionBox_RefScanMode == 2 )  bSegmentedRefLines = true; // always segmented for FLEET


  #endif

  // ---------------------------------------------------------------------------
  // calculate EPI reordering
  // ---------------------------------------------------------------------------
  #ifdef EP2D_DIFF
    switch (pMrProt->kSpace().phasePartialFourierFactor())
    {
    case SEQ::PF_OFF:
        REOInfo.usePrivatePartialFourierFactors(1.00) ;
        break ;
    case SEQ::PF_7_8:
        REOInfo.usePrivatePartialFourierFactors(0.87) ;
        break ;
    case SEQ::PF_6_8:
        REOInfo.usePrivatePartialFourierFactors(0.75) ;
        break ;
    case SEQ::PF_5_8:
        REOInfo.usePrivatePartialFourierFactors(0.66) ;
        break ;
    case SEQ::PF_HALF:
        REOInfo.usePrivatePartialFourierFactors(0.50) ; 
        break ;
    }
  #endif

  REOInfo.setModeSingleShot();
  
  
  #ifdef SUPPORT_iPAT_a_ep2d
  {
    REOInfo.setPATSegmentedRefScans(bSegmentedRefLines);
  }
  #endif

  if(! REOInfo.reorderEPI(pMrProt,pSeqLim) )mSBBErrGotoFinish(REOInfo,"REOInfo.reorderEPI");
  if(mIsTRres) REOInfo.printData(ptModule);

  // ---------------------------------------------------------------------------
  // Set gradient performance for EPI readout
  //
  // Fast gradient mode:   set slew rate to 95% of "fast" value specified 
  //                       in MeasPerm section (CHARM 324468).
  //
  // Normal gradient mode: set slew rate to 80% of "fast" value specified 
  //                       in MeasPerm section (CHARM 333764);
  //                       this option is provided to reduce Nyquist ghosts 
  //                       under low supply voltage conditions.
  // ---------------------------------------------------------------------------
  switch ( pMrProt->gradSpec().gradModeBeforeGSWD() )
  {
  case SEQ::GRAD_FAST:
      // Note that this call also sets the default gradient performance (SBBEPIReadout)
	  EPIKernel.setMinRiseTimeScaleFactor( 1.05 );
	  break;

  case SEQ::GRAD_NORMAL:
      // Note that this call also sets the default gradient performance (SBBEPIReadout)
	  EPIKernel.setMinRiseTimeScaleFactor( 1.25 );
	  break;

  default:
      textTR("Gradient mode not supported by EPI sequence");
      return SEQU_ERROR;
  }

  // ---------------------------------------------------------------------------
  // set GSWD gradient performance for excitation SBB
  // ---------------------------------------------------------------------------
  WE_121.setGSWDGradientPerformance(pMrProt, pSeqLim);

  // ---------------------------------------------------------------------------
  // configure SBBEPIKernel
  // ---------------------------------------------------------------------------
  if (pMrProt->fastImaging().freeEchoSpacing())
  {
      EPIKernel.setUseFixedEchoSpacing(pMrProt->fastImaging().echoSpacing());
  } 
  else
  {
      EPIKernel.setUseShortestEchoSpacing();
  }

  if(! EPIKernel.setPointerToReorderInfo(&REOInfo) )mSBBErrGotoFinish(EPIKernel,"EPIKernel.setPointerToReorderInfo");

  EPIKernel.setGSWDGradientPerformance        (pMrProt, pSeqLim);
  EPIKernel.setIgnoreForbiddenEchoSpacingRange(false); // important here in fSEQPrep, do not move to fSEQInit
  EPIKernel.setEchoTrainLength                (REOInfo.getEchoTrainLength());
  EPIKernel.setCenterSegment                  (REOInfo.getKSpaceCenterSegment());
  EPIKernel.setWantedTE                       (pMrProt->te()[0]);

  #ifdef SUPPORT_iPAT_a_ep2d
  {
    if (REOInfo.isPATActive() && REOInfo.getPATAccelerationFactorPE()>1)
    {
        // To achieve that the segmented reference scans have exactly the same effective
        // TE as the imaging scans we would write the following line:
        //
        // EPIKernel.setUseEchoShifting(true, REOInfo.getPATAccelerationFactorPE(), REOInfo.getPATRefCounterInSegmentWithKSCenter());
        // 
        // This leads to a non-convex reference line parameter space, because the value 
        // REOInfo.getPATRefCounterInSegmentWithKSCenter() can change depending on the
        // number of reference lines like this:
        //
        // PATRefCounterInSegmentWithKSCenter = (PATRefLinesPE/2)%PATAccelerationFactorPE
        //
        // Therefore we currently use always the minimum possible TE for the imaging scans.
        // The maximum error for the effective TE of the extra reference lines is EchoSpacing/2.
        // 
        EPIKernel.setUseEchoShifting(true, REOInfo.getPATAccelerationFactorPE(), 0);
        EPIKernel.setExpectedMaxLinesForPEBlip(REOInfo.getPATAccelerationFactorPE());
    }
    else
    {
        EPIKernel.setUseEchoShifting(false);
        EPIKernel.setExpectedMaxLinesForPEBlip(1);
    }    
  }
  #endif

  // set the number of interleaves to 1 for single-shot sequence
  // this information is required by the EPIKernel ONLY for the function calcEffEchoSpacingAndBWPerPixelPE()
  // the parameter has NO effect on the EPI readout gradient waveform
  EPIKernel.setNumInterleaves( 1 );

  /*
  #ifdef EP2D_DIFF
  if (pMrProt->PAT().PATMode()==SEQ::PAT_MODE_NONE)
    EPIKernel.setInternalPhaseCorrection                (false);
  else
    EPIKernel.setInternalPhaseCorrection                (true);
  #endif
  */


  // benpos WIP multiecho
  // set the number of echoes, and echo delay times  from specil card - times are converted from milli- to microseconds here
  EPIKernel.setNumberOfEchoes( (double)pMrProt->contrasts() ); 
  EPIKernel.setEchoDelayTimes ( l_LongEntry_delayBetweenFirstAndSecondEcho * 1000, l_LongEntry_delayBetweenRemainingEchoes * 1000 ); 
	
 //set z-shim moment from special card - gradient moment converted to mT/m*us
  EPIKernel.setZshimMomentFixed  (l_LongEntry_FixedZshim * pMrProt->te()[0] / 1000); //
  EPIKernel.setZshimMomentAltern (l_LongEntry_AlternZshim * pMrProt->te()[0] / 1000); //

  // benpos SMS-CAIPI -->
  EPIKernel.setSmsFactor(l_LongEntry_SmsFactor);
  EPIKernel.setSmsShiftFactor(l_LongEntry_SmsShiftFactor);
  // <-- benpos SMS-CAIPI


  // ggg: Setting spiral parameters
  gSpiral_setCAIPISeparationInmm(paramDoubleEchoShift);
	gSpiral_setCAIPIPeriod_us(paramDoubleCAIPIPeriod);
	gSpiral_setCAIPIDelay_us(paramDoubleCAIPIDelay);
	gSpiral_setParams(paramLongROSamples,paramLongSpiralType,
  		paramLongSpGradAmp,paramLongSpSlewRate,paramLongInterleaves,0,
  		paramLongSpiralBW,paramDoubleVD,paramDoubleAccR,paramLongZCToFlat,nIceProgram);
  	
  // ---------------------------------------------------------------------------
  // prepare SBBEPIKernel
  // ---------------------------------------------------------------------------
  if(! EPIKernel.prep(pMrProt, pSeqLim, pSeqExpo) )mSBBErrGotoFinish(EPIKernel,"EPIKernel.prep");

  // ---------------------------------------------------------------------------
  // check gradients of the EPI-kernel for grad spec violations
  // --------------------------------------------------------------------------- 
  if(! EPIKernel.checkGradients(pMrProt,pSeqLim) )
  {
      lStatus = EPIKernel.getNLSStatus();
      if (!pSeqLim->isContextPrepForBinarySearch())
          TRACE_PUT1_NLS(TC_INFO, TF_SEQ, "%s: EPIKernel.checkGradients failed with ",ptModule,EPIKernel.getNLSStatus());
      goto FINISHED;
  }

  // ---------------------------------------------------------------------------
  // check TE
  // ---------------------------------------------------------------------------
  if (EPIKernel.getNeededTE() != pMrProt->te()[0])
  {
      long   lNeededTENotOnInc = EPIKernel.getNeededTE();
      double dInc              = pSeqLim->getTE()[0].getInc();

      // ///////////////////////////////////////////////////////////////////////////////////////////////
      // #pragma message ("..........NOTE!!!!! we assume something about TE-limit-handling in MrUILink !")
      // ///////////////////////////////////////////////////////////////////////////////////////////////
      if (lNeededTENotOnInc >   10000) dInc *= 10.0;
      if (lNeededTENotOnInc > 1000000) dInc *= 10.0;

      lNeededTE = (long)(0.5+dInc*ceil(lNeededTENotOnInc/dInc));

      if ( !EPIKernel.increaseTE(lNeededTE) )mSBBErrGotoFinish(EPIKernel,"EPIKernel.increaseTE");
  }
  else
  {
      lNeededTE = pMrProt->te()[0];
  }

  // ---------------------------------------------------------------------------
  // handle setting of thickness of IR-pulse
  // ---------------------------------------------------------------------------

  // NOTE: mySeqLoop.setScaleIRThickness(2.0)) was already passed to mySeqLopp in
  // fSEQInit.
  #if (defined EP2D_DIFF) || (defined EP2D_SE)
    if(pMrProt->preparationPulses().inversion() == SEQ::SLICE_SELECTIVE)
    {
        if (! pMrProt->sliceSeries().isInterleaved())
        {
            if (!pSeqLim->isContextPrepForBinarySearch()) 
                textTR("slice selective IR only possible with interleaved excitation order");
            lStatus=SEQU_ERROR;
            goto FINISHED;
	    // ///////////////////////////////////////////////////////////////////////////////////////////////
	    // #pragma message ("..........TODO????? solve handler would be nice for conflict: IRsel <-> not interleaved excitation")
	    // ///////////////////////////////////////////////////////////////////////////////////////////////
        }
        // Reduced inversion thickness for short inversion times (e.g. STIR) - avoid crosstalk
        // for nested inversion schemes.
        if (pMrProt->ti()[0] < 500000)  // [us]
        {
            mySeqLoop.setScaleIRThickness   (1.25);
        }
    }
  #endif
  
  //benpos
  mySeqLoop.setCSatFlipAngle (l_LongEntry_FatSatFlipAngle) ;

  // ---------------------------------------------------------------------------
  // set loop parameters different from standard settings
  // --------------------------------------------------------------------------- 
  mySeqLoop.setdDistFacMinIfConcNo  (pSeqLim->getSliceDistanceFactor().getMin());
  mySeqLoop.setLinesToMeasure       (REOInfo.getLinesPerSegment()              );

  #ifdef EP2D_DIFF
    mySeqLoop.setDiffusionLoopLength(Diff->getTotalScans(pMrProt)) ;  
  #endif

  //----------------------------------------------------------------------------
  // determine number of prep-scans required
  //----------------------------------------------------------------------------
  //
  // compare SeqLoop.cpp:
  //
  // - PreparingScans = (m_lPreparingTime + (TR*Phases) - 1) / (TR*Phases)
  // - total prep scans per concatenations = PreparingScans * lSlicesInConc * Phases
  //
  //----------------------------------------------------------------------------

  if ( (pMrProt->averages() > 1) || ( pMrProt->repetitions() > 0 ) )
  {
    // we want to prepare at least for 3 sec; i.e.

    lRequiredPrepScans = (long) (3000000.0 / ( pMrProt->tr()[0]*pMrProt->physiology().phases() ) + 1.0 );
  }

  else
  {
    // no preparing for 'real' single shot
    lRequiredPrepScans = 0;
  }


  #ifdef EP2D_DIFF 

	// min. one prep scan for diffusion
	if (!lRequiredPrepScans) lRequiredPrepScans = 1;
	mySeqLoop.setPhaseCorScans(false);

  #endif
  
  #if (defined SUPPORT_iPAT_a_ep2d) && (defined PASL)

    // with iPAT a separate M0 scan BEFORE the iPAT scan is not possible
    if (REOInfo.isPATActive() && REOInfo.getPATAccelerationFactorPE()>1)
    {
        bEnableFirstPrepScanAsM0Scan = false;
    }

  #endif

  #ifdef PASL 

    if ( pMrProt->repetitions() > 0 )
    {
      // minimum PrepScanTime [ms] is taken from PaslUI
      lRequiredPrepScans = (long)( ceil(Pasl.getPrepScanTime()*1e3 / pMrProt->tr()[0]) );
    }
    else
    {
      // no preparing for 'real' single shot
      lRequiredPrepScans = 0;

      // no separate M0 scan for single shot
      bEnableFirstPrepScanAsM0Scan = false;
    }

    if( bEnableFirstPrepScanAsM0Scan )
    {
      // ensure at least a single M0 scan: >= 1;
      lRequiredPrepScans = max(1, lRequiredPrepScans);

      // decrease repetitions by one (first repetition is M0 reference scan measured as PrepScan)
      mySeqLoop.setRepetitionsToMeasure( pMrProt->repetitions() - 1 );
    }

  #endif

  #ifdef SUPPORT_iPAT_a_ep2d

	m_lInitialDummyScans = 0;
	m_lPostDummyScans    = 0;
	m_lSliceAccPrepScans = 0;
	
	if (REOInfo.isPATActive() && REOInfo.getPATAccelerationFactorPE()>1)
	{
		//benpos 
		if (l_SelectionBox_RefScanMode == 0 || l_SelectionBox_RefScanMode == 1 ) // normal EPI like reference scans (segmented or single-shot)
		{
		   
		   m_lInitialDummyScans = lRequiredPrepScans;
		   if (l_LongEntry_SmsFactor > 1.0)
			{
				m_lPostDummyScans    = lRequiredPrepScans;
			}
			else
			{
				m_lPostDummyScans    = 0;				
				#ifdef PASL
				m_lPostDummyScans = lRequiredPrepScans; 
				#endif
			}

						
			if (lMinPrepScansNoPATRefScans>m_lInitialDummyScans)
			{
				m_lInitialDummyScans = lMinPrepScansNoPATRefScans;
				lRequiredPrepScans = lMinPrepScansNoPATRefScans;
			}
			if (bSegmentedRefLines)
			{
				m_lPATPrepScans =  REOInfo.getPATAccelerationFactorPE();
			}
			else
			{
				m_lPATPrepScans = 1;
			}
				
			
		
		}
		else if (l_SelectionBox_RefScanMode == 2  ) // FLEET refscans
		{
		   
		   m_lInitialDummyScans = 0;
		   m_lPostDummyScans    = lRequiredPrepScans;

		   m_lPATPrepScans =  REOInfo.getPATAccelerationFactorPE() + l_LongEntry_FleetDummies;
	
		}
		else // SBBPatRefScan like reference lines
		{
			
			m_lInitialDummyScans = 0;
			m_lPostDummyScans    = lRequiredPrepScans; 

			//lPrepScansNoPATRefScans = lRequiredPrepScans; 
			SBBPATRefScan.setlMatchedCenterLine      (REOInfo.getKSCenterLin());  // m_lMatchedCenterLine/PartForPATRefScan may be set by sequence (default is -1, i.e. sequence didn't set anything -> SBBPATRefScan will not match center line)
			SBBPATRefScan.setlMatchedCenterPartition (REOInfo.getKSCenterPar());  // MatchedCenterLine/Par: if set, the line/part. numbering of the refscan will get an offset so that refscan has the same centerLin/Par. number as 'main' scan (depends on the ICE program if this is required)
																					// 	m_lMatchedCenterPartForPATRefScan	
			
	//		SBBPATRefScan.setlBaseRes(long(d_DoubleEntry_SSBPatRefscanBaseRes) ) ;  // pMrProt->kSpace().baseResolution());
	//		SBBPATRefScan.setlBandwidth (long(d_DoubleEntry_SSBPatRefscanBW));
	//		SBBPATRefScan.setdTE (d_DoubleEntry_SSBPatRefscanTE);	
			SBBPATRefScan.setbUseMinTE(true); 
			if( !SBBPATRefScan.prep(pMrProt, pSeqLim, pSeqExpo)) 
			 {
			   if ( ! pSeqLim->isContextPrepForBinarySearch() )
			   {
				   TRACE_PUT1(TC_ALWAYS, TF_SEQ,"%s: Error encountered in preparing SBBPATRefScan", ptModule );
			   }
			   lStatus = SBBPATRefScan.getNLSStatus();
			   goto FINISHED;
			}
			SBBPATRefScan.setRequestsPerMeasurement( SBBPATRefScan.getlPartitionsToMeasure() * pMrProt->sliceSeries().size() );
		
		} //end what type refscan mode 
	
	} 
	else // iPAT is off:
	{
		if (l_LongEntry_SmsFactor > 1.0)
		{
			m_lInitialDummyScans =0;
			m_lPostDummyScans    = lRequiredPrepScans;
		}
		else
		{
			m_lInitialDummyScans = lRequiredPrepScans;
			m_lPostDummyScans    = 0;

		}

	}

	
	if (l_LongEntry_SmsFactor > 1.0)
	{
		m_lSliceAccPrepScans = 1;
	}
	else
	{
		m_lSliceAccPrepScans = 0; 
	}


  #endif
  
  
// benpos SMS-CAIPI -->

//calculate values for the prescans 
  
  #ifdef TGSE
   lRequiredPrepScans = m_lInitialDummyScans;
   EPIKernel.set_CycleLength(2 * rAsl.getulArrayLength());
  #else
  // Add up total number of prep scans
  lRequiredPrepScans = m_lInitialDummyScans; //remove for FLEET - ew just keep adding to th eoriginal value! 
  #endif

#ifdef SUPPORT_iPAT_a_ep2d
  lRequiredPrepScans += m_lPATPrepScans;
#endif

#ifdef COMPILE_EP2D_DIFF
  lRequiredPrepScans += m_lAdjPrepScans;
#endif

  
  lRequiredPrepScans += m_lPostDummyScans;

  lRequiredPrepScans += m_lSliceAccPrepScans;
 
  // ORDER of prep scans
  //
  // m_lInitialDummyScans-----m_lSliceAccPrepScans-----m_lPATPrepScans-------m_lInitialDummyScans-------m_lAdjPrepScans--------Imaging-------...
  //
  // ----------1----------------------1---------------------1/0----------------------1---------------------DAQ-----------------DAQ---------...
  // 
  // --------dummy--------------------DAQ-------------------DAQ---------------------dummy------------------DAQ-----------------DAQ---------...
  // 
  // --------all slices------------all slices------------all slices------------reduced slices---------reduced slices------reduced slices---...

 // <-- benpos SMS-CAIPI


  // ---------------------------------------------------------------------------
  // force SeqLoop to use a specific time for prep scans
  // ---------------------------------------------------------------------------
  mySeqLoop.setlPreparingTime(lRequiredPrepScans*pMrProt->tr()[0]*pMrProt->physiology().phases());

  // ---------------------------------------------------------------------------
  // Advise SeqLoop to perform prep-scans only within the first measurement.
  // ---------------------------------------------------------------------------
  if (bPrepScansOnlyInFirstMeasurement)
  {
      mySeqLoop.setePerformPreparingScans(OnlyFirstRepetition);
  }
  else
  {
      mySeqLoop.setePerformPreparingScans(Always);
  }

  // ---------------------------------------------------------------------------
  // Set long TR triggering mode
  // ---------------------------------------------------------------------------

  // switch on by default

  mySeqLoop.setLongTRTrigMode( true );

  // no long TR triggering mode if navigator triggering is switched on

  #ifdef SUPPORT_PACE
  
	if ( pMrProt->NavigatorParam().RespComp() != SEQ::RESP_COMP_OFF )
	{
		mySeqLoop.setLongTRTrigMode( false );
	}
  
	else
	{    
		mySeqLoop.setLongTRTrigMode( true );
	}
  
  #endif
	  
  // no long TR triggering mode for ep2d_se, ep2d_fid and ep2d_pasl sequences
	  
  #if ( (defined EP2D_SE) || (defined EP2D_FID) || (defined PASL) )

	mySeqLoop.setLongTRTrigMode( false );

  #endif

  // ---------------------------------------------------------------------------
  // CHARM 305893:
  // 
  // TR should ALWAYS be the physically correct TR for single shot EPI. Also and
  // especially over multiple measurements. This leads to the following restrictions:
  //
  // - multi-slice mode must be interleaved
  //
  // - We do not allow repetition delay times, because this would disturb the
  //   steady state and change the true TR. This restriction can be overcome
  //   with a new parameter which allows to force SeqLoop to use a certain TR-fill
  //   at the end of the concatenation. So the measurement delay is included in
  //   TR.
  //
  // - For the same reason the measPause must be zero.
  //
  // - It is not allowed to use multiple concatenations and multiple measurements
  //   at the same time. Again this would lead to physically incorrect TRs for the
  //   first acquisition of a certain slice within the measurement for all measurements
  //   except the first, i.e. for all repetitions.
  // ---------------------------------------------------------------------------
  if (pMrProt->kSpace().multiSliceMode() != SEQ::MSM_INTERLEAVED)
  {
      if (!pSeqLim->isContextPrepForBinarySearch()) textTR("pMrProt->kSpace().multiSliceMode() == SEQ::MSM_INTERLEAVED required!");
      lStatus=SEQU_ERROR;
      goto FINISHED;
  }

  for (lRepetition = 0; lRepetition < K_NO_REP_TIMES; lRepetition++)
  {   
      if (pMrProt->repetitionDelayTime()[lRepetition])
      { 
          if (!pSeqLim->isContextPrepForBinarySearch()) textTR("All pMrProt->repetitionDelayTime()[x] must be zero!");
          lStatus=SEQU_ERROR;
          goto FINISHED;
      }
  }

  if (pMrProt->measPause())
  {
      if (!pSeqLim->isContextPrepForBinarySearch()) textTR("pMrProt->measPause() must be zero!");
      lStatus=SEQU_ERROR;
      goto FINISHED;
  }

  // -------------------------------------------------------------------------
  // Apply restrictions for multiple concatenations
  // -------------------------------------------------------------------------
  if ( pMrProt->concatenations() > 1 )
  {		
	  bool bMultiConcatsAllowed = false;
	
	  // Multiple concatenations are allowed if long TR triggering mode is enabled
	  // and standard triggering is active

	  if ( mySeqLoop.isLongTRTrigMode() )
	  {
		  SEQ::PhysioSignal FirstSignal;
		  SEQ::PhysioMethod FirstMethod;
		  SEQ::PhysioSignal SecondSignal; 
		  SEQ::PhysioMethod SecondMethod;

		  pMrProt->physiology().getPhysioMode (FirstSignal, FirstMethod, SecondSignal, SecondMethod);
	  
		  if ( FirstMethod == SEQ::METHOD_TRIGGERING )
		  {		
			  bMultiConcatsAllowed = true;
		  }
	  }
			
	  // Multiple concatenations are allowed if navigator triggering is active

	  #ifdef SUPPORT_PACE

		  if ( pMrProt->NavigatorParam().RespComp() != SEQ::RESP_COMP_OFF ) 
			  bMultiConcatsAllowed = true;

	  #endif

	  // Multiple concatenations are allowed with ep2d_se and ep2d_fid 
	  // sequences if there is a single repetition

	  #if ( (defined EP2D_SE) || (defined EP2D_FID) )

		  if ( pMrProt->repetitions() == 0 ) 
			  bMultiConcatsAllowed = true;

	  #endif

	  // Return with error if multiple concatenations are not allowed

	  if ( !bMultiConcatsAllowed )
	  {	
		  lStatus = SEQU_ERROR;
	  	
		  if (! pSeqLim->isContextPrepForBinarySearch() )  	
		  {	  	
			  TRACE_PUT1(TC_INFO, TF_SEQ, "%s: Multiple concatenations are not possible ", ptModule );	
		  }

		  goto FINISHED;
	  }
  }

  // navigator triggering

    if ((pMrProt->NavigatorParam().RespComp() == SEQ::RESP_COMP_TRIGGER)
        || (pMrProt->NavigatorParam().RespComp()  == SEQ::RESP_COMP_TRIGGER_AND_FOLLOW)) 
    {
        mySeqLoop.setOptfsPrepflag(true);
    }
    else
    {
        mySeqLoop.setOptfsPrepflag(false);
    }

  // ---------------------------------------------------------------------------
  // prepare standard loop     
  // ---------------------------------------------------------------------------
  if(! mySeqLoop.prep(pMrProt, pSeqLim, pSeqExpo) )mSBBErrGotoFinish(mySeqLoop,"mySeqLoop.prep");

  #ifndef PASL    
    #ifdef PACE3D
    mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(TWAKEUP);
    #else
    mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(0);
    #endif
  #else

  if (pMrProt->preparationPulses().fatSatMode() == SEQ::FAT_SAT_STRONG)
  {
      double dRfEnergyInSBBs = 0.0;

	  
	  if ( pMrProt->repetitions() > 0 )
      {
        // FatSat / Spoil gradient preparation - FATSAT_ALL_SLICES
        CSatFat.setRequestsPerMeasurement  ( mySeqLoop.getKernelRequestsPerMeasurement(pSeqLim,SecondMeas) );
      }
      else
      {
        // calc request for FirstMeas without PrepScans
        CSatFat.setRequestsPerMeasurement  ( mySeqLoop.getKernelRequestsPerMeasurement(pSeqLim,FirstMeas) /    (lRequiredPrepScans + 1) );
      }
      CSatFat.setCSatMode                ( SBBCSatCode_Fat  );
      CSatFat.setIdent                   ( "PaslFS2"         );
      CSatFat.setGSWDGradientPerformance ( pMrProt, pSeqLim );
	  CSatFat.setFlipAngle(l_LongEntry_FatSatFlipAngle); //benpos add this here.
	  
	  SpoilGrad.setRequestsPerMeasurement  ( CSatFat.getRequestsPerMeasurement() );
      SpoilGrad.setGSWDGradientPerformance ( pMrProt, pSeqLim );
    
      if(!SBB.prepSBBAll(pMrProt,pSeqLim,pSeqExpo,&dRfEnergyInSBBs) )
	  { 
		  return (SBB.getpSBBLastPrep()->getNLSStatus());
	  }
  }
  // set Requests for FatSat in Pasl
  Pasl.setRequestsPerMeasurement  ( 1 );
  // set time between start of SBB and excitation of the first image slice
  Pasl.setTimeToImageSliceExc_us  ( WE_121.getlRFCenterTime() );

  //benpos set PASL fatsat flipangle (the first FatSat is played out by PASL sbb)
  Pasl.setCSatFlipAngle (l_LongEntry_FatSatFlipAngle);
    
  if( !Pasl.prep(pMrProt, pSeqLim, pSeqExpo) )
  {
	  return (Pasl.getNLSStatus()); 
  }
 
  cout << " SBB  ASL mode is " << Pasl.getASLMode() << endl;
cout << " pMrProt->Asl().Mode() " << pMrProt->Asl().Mode() <<  endl;

	   

  if (pMrProt->preparationPulses().fatSatMode() == SEQ::FAT_SAT_STRONG)
  {    
  #ifdef PACE3D
	  // FATSAT_ALL_SLICES
      // tell seqloop that first FatSat is not to be taken into account
      mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(TWAKEUP + Pasl.getDurationPerRequest() 
        - CSatFat.getDurationPerRequest() - SpoilGrad.getDurationPerRequest() );
  #else
        // FATSAT_ALL_SLICES
        // tell seqloop that first FatSat is not to be taken into account
      //mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(Pasl.getDurationPerRequest()  - CSatFat.getDurationPerRequest() - SpoilGrad.getDurationPerRequest() )  ;
	  //benpos: additionally we need to substract the time for non-acquired slices when using SMS
	  mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(Pasl.getDurationPerRequest()  - CSatFat.getDurationPerRequest() - SpoilGrad.getDurationPerRequest()  - (l_LongEntry_SmsFactor-1)/l_LongEntry_SmsFactor*pMrProt->sliceSeries().size() ) ;
  #endif
  }
  else
  {
  #ifdef PACE3D
	  mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(TWAKEUP + Pasl.getDurationPerRequest() );
  #else
	  mySeqLoop.setAdditionalTimeForExtraEBinTRUsec(Pasl.getDurationPerRequest());
  #endif
  }
  #endif


  //----------------------------------------------------------------------------
  // Check, if we got the correct number of preparing scans from SeqLoop.
  // This is important for PAT-measurements, because reference lines are acquired
  // during the last prep-scans.
  //----------------------------------------------------------------------------
  #ifdef SUPPORT_iPAT_a_ep2d
  {
    if (mySeqLoop.getlPreparingScans() != lRequiredPrepScans)
    {
        if (!pSeqLim->isContextPrepForBinarySearch()) TRACE_PUT3(TC_INFO, TF_SEQ, "%s: NEED %ld PrepScans from SeqLoop, BUT got %ld",ptModule,lRequiredPrepScans,mySeqLoop.getlPreparingScans());                   
        lStatus=SEQU_ERROR;
        goto FINISHED;
    }
  }
  #endif
  
  // ---------------------------------------------------------------------------
  // calculate the TR/TI fill times
  // --------------------------------------------------------------------------- 
  if (! calculateTRTIFillTimes(pMrProt,pSeqLim,pSeqExpo,&lNeededTI,&lNeededTR) )
  {
      if (!pSeqLim->isContextPrepForBinarySearch()) {CheckStatusPB(SEQU_ERROR, "calculateTRTIFillTimes failed");}
      else                                          {CheckStatusB (SEQU_ERROR);}
  }

  // ---------------------------------------------------------------------------
  // prep pace stuff:
  // copy slice position data
  // --------------------------------------------------------------------------- 
  #ifdef PACE3D
   if ( !gPaceFeedback.Prep(pMrProt,asSLC) )
    {
      // error can only be caused by programming error
        // => trace also if pSeqLim->isContextPrepForBinarySearch()
        textTR("PaceFeedback.Prep failed.");
        lStatus=SEQU_ERROR;
        goto FINISHED;
    }            
  #endif 

  // ---------------------------------------------------------------------------
  // calculate energy
  // ---------------------------------------------------------------------------
  dEnergy  = mySeqLoop.getEnergy(pMrProt);
  // benpos SMS-CAIPI -->
  if (l_LongEntry_SmsFactor > 1.0)
  {
      lFullTRMeas = 0;	// number of measurements at full non-sliceaccelerated TR for prep scans (anything before slice acc scans)
    #ifdef SUPPORT_iPAT_a_ep2d
          lFullTRMeas += m_lPATPrepScans;
    #endif
      lFullTRMeas += m_lSliceAccPrepScans;
      lFullTRMeas += m_lInitialDummyScans;

	  // m_lInitialDummyScans-----m_lSliceAccPrepScans-----m_lPATPrepScans-------m_lInitialDummyScans-------m_lAdjPrepScans--------Imaging-------...
 	  TRACE_PUT5(TC_INFO, TF_SEQ, "%s: m_lInitialDummyScans: %ld, m_lSliceAccPrepScans: %ld, m_lPATPrepScans: %ld, lRequiredPrepScans: %ld", ptModule, m_lInitialDummyScans, m_lSliceAccPrepScans, m_lPATPrepScans, lRequiredPrepScans);

      lTotalMeas = mySeqLoop.getKernelRequestsPerMeasurement(pSeqLim,FirstMeas ) / pMrProt->sliceSeries().size();	// total number of measurements (TRs) for prep scans and slice acc scans (excluding repetitions)

      lReducedTRMeas = lTotalMeas - lFullTRMeas;	// number of measurements (TRs) for slice accelerated scans 
      lReducedTRMeas_SecondMeas = (mySeqLoop.getKernelRequestsPerMeasurement(pSeqLim,SecondMeas) / pMrProt->sliceSeries().size()) * pMrProt->repetitions();	// number of repetitions

      dEnergySeqLoopPerSlice = dEnergy/(double)((lTotalMeas+lReducedTRMeas_SecondMeas)* pMrProt->sliceSeries().size());

      long lNSlicesAcc = pMrProt->sliceSeries().size() / l_LongEntry_SmsFactor;
      dEnergy = 0.0;
      dEnergy += dEnergySeqLoopPerSlice * lReducedTRMeas * lNSlicesAcc;
      dEnergy += dEnergySeqLoopPerSlice * lFullTRMeas * pMrProt->sliceSeries().size();
      
      //energy of full TR scans
      EPIKernel.setRequestsPerMeasurement(lFullTRMeas * pMrProt->sliceSeries().size());
      dEnergy += EPIKernel.getEnergyPerMeasurement(pMrProt->sliceSeries().size());

	  //energy in slice accelerated scans
      EPIKernel.setRequestsPerMeasurement(lReducedTRMeas * lNSlicesAcc);
      dEnergy += EPIKernel.getEnergyPerMeasurementSms(lNSlicesAcc);

      EPIKernel.setRequestsPerMeasurement(lReducedTRMeas_SecondMeas * lNSlicesAcc);
      dEnergy += EPIKernel.getEnergyPerMeasurementSms(lNSlicesAcc);
      
	  dEnergy += dEnergySeqLoopPerSlice * lReducedTRMeas_SecondMeas * lNSlicesAcc;
  }
  else // <-- benpos SMS-CAIPI  
  {    // this is the regular energy calculation
	  EPIKernel.setRequestsPerMeasurement(mySeqLoop.getKernelRequestsPerMeasurement(pSeqLim,FirstMeas ));
	  dEnergy += EPIKernel.getEnergyPerMeasurement(pMrProt->sliceSeries().size());
	  EPIKernel.setRequestsPerMeasurement(mySeqLoop.getKernelRequestsPerMeasurement(pSeqLim,SecondMeas));
	  dEnergy += EPIKernel.getEnergyPerMeasurement(pMrProt->sliceSeries().size()) * pMrProt->repetitions();
  }

  
  
  //benpos sbbpatrefscan: add energy for patrefscan
  if (l_SelectionBox_RefScanMode==3)
  {
	  cout << " TODO: benpos, you need to fix the prepscan counter and energy calculations for SBBParRefScan! " << endl;
	  dEnergy+= SBBPATRefScan.getEnergyPerRequest() * pMrProt->sliceSeries().size();
  }


  #ifdef PASL
  
  dEnergy += ( Pasl.getRequestsPerMeasurement() * Pasl.getEnergyPerRequest() )   * ( pMrProt->repetitions() + 1 + lRequiredPrepScans )  / l_LongEntry_VolumesPerInversion; //benpos

  //Add energy of Strong FatSat
  if( pMrProt->preparationPulses().fatSatMode() == SEQ::FAT_SAT_STRONG )
  {
    if( (CSatFat.getRequestsPerMeasurement() > 1) )
    {
      dEnergy += ( CSatFat.getRequestsPerMeasurement() - 1 ) * CSatFat.getEnergyPerRequest() * ( pMrProt->repetitions() + 1 + lRequiredPrepScans )   / l_LongEntry_VolumesPerInversion; //benpos
    }
  }

  if( bEnableFirstPrepScanAsM0Scan )
  {
    // one volume is explicitly taken out of SeqLoop (pMrProt->repetitions() - 1) and acquired as PrepScan (M0 volume)
    EPIKernel.setRequestsPerMeasurement(mySeqLoop.getKernelRequestsPerMeasurement(pSeqLim,SecondMeas));
    dEnergy -= EPIKernel.getEnergyPerMeasurement(pMrProt->sliceSeries().size());
    
    // subtract one volume in SeqLoop and one REFERENCE (RF switched off) from PrepScans
    dEnergy -= 2 * Pasl.getRequestsPerMeasurement() * Pasl.getEnergyPerRequest();

    // subtract energy of Strong FatSat (taken out of SeqLoop)
    if (pMrProt->preparationPulses().fatSatMode() == SEQ::FAT_SAT_STRONG)
    {
      if ( CSatFat.getRequestsPerMeasurement() > 1)
      {
        dEnergy -= ( CSatFat.getRequestsPerMeasurement() - 1 ) * CSatFat.getEnergyPerRequest();
      }
    }
  }

  #endif 
  //---------------------------------------------------------------------------
  // activate online ice-process
  //---------------------------------------------------------------------------
  // ggg
  if(SpiralModeOn()) {
  	EPIKernel.getReadOutAddress()->Mdh.deleteFromEvalInfoMask(MDH_ONLINE);
  	EPIKernel.getReadOutAddress()->Mdh.addToEvalInfoMask(MDH_OFFLINE);
  } else {
  	EPIKernel.getReadOutAddress()->Mdh.addToEvalInfoMask(MDH_ONLINE);
  }
  

  // ---------------------------------------------------------------------------
  // set the gain of the receiver
  // ---------------------------------------------------------------------------
  lStatus = fSSLSetRxGain(K_RX_GAIN_CODE_HIGH,pMrProt,pSeqLim);

  // error can only be caused by programming error
  // => trace also if pSeqLim->isContextPrepForBinarySearch()
  CheckStatusPB(lStatus,"fSSLSetRxGain");

  // ---------------------------------------------------------------------------
  // calculate effective echo-spacing 
  // and bandwidth per pixel in phase-encode direction
  // ---------------------------------------------------------------------------
  if (! EPIKernel.calcEffEchoSpacingAndBWPerPixelPE( pMrProt, lEffectiveEchoSpacing, dBandwidthPerPixelPE ) )
  {
	  TRACE_PUT1(TC_INFO, TF_SEQ, "%s ERROR: calcEffEchoSpacingAndBWPerPixelPE() failed", ptModule );            	
	  lStatus=SEQU_ERROR;
	  goto FINISHED;
  }



// benpos calculate measurement time

if ( l_LongEntry_SmsFactor == 1 &&  REOInfo.getPATAccelerationFactorPE() == 1 )   
{
	 
	pSeqExpo->setTotalMeasureTimeUsec  (      pMrProt->tr()[0] / l_LongEntry_VolumesPerInversion * (pMrProt->repetitions()+1) 
										  +  m_lInitialDummyScans * pMrProt->tr()[0] // here InitialDummyScans have the full TR, no other prepscans
										);
}
else 
{
	pSeqExpo->setTotalMeasureTimeUsec  (      pMrProt->tr()[0] / l_LongEntry_VolumesPerInversion * (pMrProt->repetitions()+1) 
										  +  m_lInitialDummyScans * (mySeqLoop.getlKernelScanTime() / l_LongEntry_SmsFactor   ) 
								//dimo	  +  m_lSliceAccPrepScans * (mySeqLoop.getlKernelScanTime()
										  +  m_lPATPrepScans * (mySeqLoop.getlKernelScanTime() /  l_LongEntry_SmsFactor ) 
										  +  m_lPostDummyScans * pMrProt->tr()[0] // only PostDummyScans have the full TR	
										);
}


  // ---------------------------------------------------------------------------
  // prepare common exports
  // ---------------------------------------------------------------------------
  pSeqExpo->clearScanningSequence           ();
  pSeqExpo->addtoScanningSequence           ("EP");
  pSeqExpo->setPhaseCorScans                (2                                                  );
  //pSeqExpo->setTotalMeasureTimeUsec         ( (mySeqLoop.getTotalMeasTimeUsec  (pMrProt, pSeqLim) - (m_lInitialDummyScans + m_lPATPrepScans + m_lPostDummyScans + m_lSliceAccPrepScans) * Pasl.getDurationPerRequest()) / l_LongEntry_VolumesPerInversion + (m_lInitialDummyScans + m_lPATPrepScans + m_lPostDummyScans + m_lSliceAccPrepScans) * mySeqLoop.getlKernelScanTime() );
  //pSeqExpo->setTotalMeasureTimeUsec         ( pMrProt->tr()[0] / l_LongEntry_VolumesPerInversion * (pMrProt->repetitions()+1) );
  pSeqExpo->setRFEnergyInSequence_Ws        (dEnergy                                            );
  pSeqExpo->setMeasuredPELines              (REOInfo.getLinesToMeasure()                        );
  pSeqExpo->setMeasured3dPartitions         (REOInfo.getPartitionsToMeasure()                   );
  pSeqExpo->setEchoSpacing                  (EPIKernel.getEchoSpacing()                         );
  pSeqExpo->setEffectiveEpiEchoSpacing      (lEffectiveEchoSpacing                              );
  pSeqExpo->setBandwidthPerPixelPhaseEncode (float(dBandwidthPerPixelPE)                        );
  pSeqExpo->setPCAlgorithm                  (SEQ::PC_ALGORITHM_NONE                             );
  pSeqExpo->setRelevantReadoutsForMeasTime  (mySeqLoop.getNumberOfRelevantADCs()                );
  pSeqExpo->setB0Correction                 (pMrProt->NavigatorParam().RespComp() == SEQ::RESP_COMP_OFF);
  pSeqExpo->setMaxwellCorrection            (true                                               );
  pSeqExpo->setMaxwellIntegralROGradient    (float(EPIKernel.getGRO().getMaxwellIntegral())     );

  // ---------------------------------------------------------------------------
  // prepare exports concerning measurement time per mesurement
  // ---------------------------------------------------------------------------
  if (pMrProt->repetitions())
  {
      pSeqExpo->setMeasureTimeUsec             ( mySeqLoop.getMeasurementTimeUsec(pMrProt, pSeqLim,SecondMeas));
      pSeqExpo->setPreparingTimeInFirstMeasUSec( mySeqLoop.getMeasurementTimeUsec(pMrProt, pSeqLim,FirstMeas)
                                                -mySeqLoop.getMeasurementTimeUsec(pMrProt, pSeqLim,SecondMeas)
                                               );
  }
  else
  {
      pSeqExpo->setMeasureTimeUsec             (mySeqLoop.getMeasurementTimeUsec(pMrProt, pSeqLim,FirstMeas));
      pSeqExpo->setPreparingTimeInFirstMeasUSec(0);
  }

  // ---------------------------------------------------------------------------
  // determine default ice-program
  // ---------------------------------------------------------------------------
  bSuccess = mySeqLoop.setIceProgram 
             (
                pMrProt, pSeqLim, pSeqExpo,
                mySeqLoop.getDataRateMBytePerSec(pMrProt, pSeqLim, REOInfo.getEchoTrainLength()),
                SEQ::RS_LIN_IN_PAR
             );
  if (!bSuccess) mSBBErrGotoFinish(mySeqLoop,"mySeqLoop.setIceProgram failed");

  // ---------------------------------------------------------------------------
  // activate zero-order, first-order and primary mode 
  // phase-correction in ICE program
  // ---------------------------------------------------------------------------

  #ifndef EP2D_DIFF
  pSeqExpo->setOnlinePhaseCorrectionAlgo(ICE_ONLINEPC_AUTOCORR|ICE_ONLINEPC_CROSSCORR|ICE_ONLINEPC_PRIMARYMODE );
  #else
  pSeqExpo->setOnlinePhaseCorrectionAlgo(ICE_ONLINEPC_AUTOCORR|ICE_ONLINEPC_CROSSCORR_ACROSSSEGMENTS_EPI|ICE_ONLINEPC_PRIMARYMODE);
  #endif


  // benpos WIP - overwrite the setting for phase correction algorythm for using the new "local" correction if selected if this is selected in protocol (l_SelectionBox_PhaseCorrMode = 0, the default) 
  // this uses a different PC approach in OnlineTSE
  // this modification is chosed after discusion with Thorsten Feiweier @ Siemens Erlangen
  // it should reduce the increasing ghosting over time when using the 32ch coils in long scans
 
  if(l_SelectionBox_PhaseCorrMode==1)  pSeqExpo->setOnlinePhaseCorrectionAlgo(ICE_ONLINEPC_AUTOCROSSCORR_ACROSSSEGMENTS | ICE_ONLINEPC_PRIMARYMODE);

  // to collapse back to pixel-by-pixel phase conjugation
  if(l_SelectionBox_PhaseCorrMode==2)  pSeqExpo->setOnlinePhaseCorrectionAlgo(0 );



  // ---------------------------------------------------------------------------
  // variant specific exports, ice program selection, etc.
  // ---------------------------------------------------------------------------
  #ifdef EP2D_DIFF
  pSeqExpo->setAdditionalScaleFactor(1);
  pSeqExpo->setICEProgramParam(ICE_PROGRAM_PARA_ICE_MEM_MB,0);

    if (pMrProt->diffusion().mode()==SEQ::DIFFMODE_TENSOR ||
        (pMrProt->diffusion().mode()==SEQ::DIFFMODE_FREE && pMrProt->diffusion().diffDirections() >= 6) ) 
    {
        if (pMrProt->NavigatorParam().RespComp() != SEQ::RESP_COMP_OFF )
        {
            if( pSeqLim->isContextNormal() ) TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d)",__FILE__,__LINE__);
            return SEQU_ERROR;
        }

        if( pMrProt->sliceGroupList().size() > 1 ) // if DTI active -> Mosaic could have been selected. MultiSliceMultiAngle is not allowed then.
        {
            if( pSeqLim->isContextNormal() ) TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d)",__FILE__,__LINE__);
            return SEQU_ERROR;
        }
      #ifdef WIP
        pSeqExpo->setICEProgramFilename("%CustomerIceProgs%\\IceProgramDti2D") ;
      #else
        pSeqExpo->setICEProgramFilename("%SiemensIceProgs%\\IceProgramDti2D") ;
      #endif
	pSeqExpo->setNSet               (Diff->getTotalScans(pMrProt)*pMrProt->measurements()) ;
	pSeqExpo->setOnlineFFT          (SEQ::ONLINE_FFT_PHASE);
	pSeqExpo->setICEProgramParam    (ICE_PROGRAM_PARA_SHOW_OFFLINE, SEQ::SO_SHOW_NO);

	// set flag to exclude SET and REP dimensions from data amount calculation in Sequence.cpp
	// this is appropriate for DTI scans, which perform inline image calculation without raw data storage
	pSeqLim->setInteractiveRealtime (SEQ::ON);
    }
    else {
      #ifdef WIP
        pSeqExpo->setICEProgramFilename("%CustomerIceProgs%\\IceProgramDiffusion2D");
      #else
        if( pMrProt->NavigatorParam().RespComp() != SEQ::RESP_COMP_OFF && pMrProt->NavigatorParam().RespComp() != SEQ::RESP_COMP_BREATH_HOLD )
        {
            pSeqExpo->setICEProgramFilename("%SiemensIceProgs%\\IceProgramDiffusion2D+PACE");
        }
        else
        {
            pSeqExpo->setICEProgramFilename("%SiemensIceProgs%\\IceProgramDiffusion2D");
        }
        if( pMrProt->averages() > 1 )
        {
            //  Overwrite setting from SBBDiffusion_Bipolar based on # repetitions
            pSeqExpo->setICEProgramParam(ICE_PROGRAM_DIFF_SCALING,pMrProt->averages());
            //  The extract functor devides AdditionalScaleFactor by # averges
            //  Scale due to multiple averages is already considered by diffusion functor
            pSeqExpo->setAdditionalScaleFactor(pMrProt->averages());

            //  Calculate size of raw data object (see Sequence::checkRawDataSizeOfMeasurement)
            
            double dMemoryUsedByImageRecon_MByte = (double) pMrProt->averages()
                * (double)pMrProt->sliceSeries().size()
                * (double)maximum( 1L, pMrProt->physiology().phases() )
                * (double)pMrProt->contrasts() //Ecos
                * (double)maximum( 1L, pMrProt->combinedEchoes() )
                * (double)(pMrProt->repetitions()+1)
                * (double)pMrProt->kSpace().baseResolution()
                * (double)pMrProt->coilInfo().Meas().getNumOfUsedRxChan()
                * (double)Diff->getTotalScans(pMrProt)
                * 8.0 / 1024.0 / 1024.0;
            
     
            if(pMrProt->kSpace().Interpolation2D())
            {
                // Only *2, because ROFT is done in raw data object. In PE direction no FT is done in the raw data object and
                // has not to be taken into account.
                dMemoryUsedByImageRecon_MByte *= 2;
            }
            
            //---------------------------------------------------------------------------
            // Calculate the additional memory demand considering lines and partitions
            //---------------------------------------------------------------------------
            int i3DFTLength = 1;
            if( pSeqExpo->getMeasured3dPartitions() > 1 ) //NoOfFourierPartitions > 1 == 3D case
            {
                i3DFTLength = (long)(pMrProt->kSpace().imagesPerSlab() * (pMrProt->kSpace().sliceOversampling() + 1.0));
                
                if( pMrProt->kSpace().sliceResolution() > 1.0 )
                    i3DFTLength *= int(pMrProt->kSpace().sliceResolution());
                
                i3DFTLength = NINT(i3DFTLength / 2.0) * 2;
            }
             
            const SEQ::PATSelMode PATMode = pMrProt->PAT().PATMode();
            int iNReflin = 0;
            int iAccelFact = 1;
            //---------------------------------------------------------------------------
            // Take the PAT mode into account
            //---------------------------------------------------------------------------
            if(PATMode != SEQ::PAT_MODE_NONE)
            {
                iNReflin = pMrProt->PAT().lRefLinesPE + 1;  // regard also iPAT reference lines
                // + 1: to be sure that rawdata is not to much because division below may truncate to next smaller integer.
                iAccelFact = pMrProt->PAT().AccelFactPE() * pMrProt->PAT().AccelFact3D();
            }
            dMemoryUsedByImageRecon_MByte *= ((double)(pSeqExpo->getMeasuredPELines() / iAccelFact + iNReflin) * (double)i3DFTLength);
            pSeqExpo->setICEProgramParam(ICE_PROGRAM_PARA_ICE_MEM_MB,long(dMemoryUsedByImageRecon_MByte+0.5));
        }
      #endif

        pSeqExpo->setNSet               (Diff->getTotalScans(pMrProt));
	// set flag to include SET and REP dimensions in data amount calculation in Sequence.cpp
	// this is appropriate for non-DTI scans, which save raw data and do not perform inline image calculation
	pSeqLim->setInteractiveRealtime (SEQ::OFF);
    }

    pSeqExpo->setPCAlgorithm         (SEQ::PC_ALGORITHM_NONE);    

    pSeqExpo->setApplicationCard     (SEQ::APPLICATION_CARD_DIFF     );
    pSeqExpo->setApplicationCardName (SEQ::APPLICATION_CARD_NAME_DIFF);

    if(pMrProt->preparationPulses().inversion() == SEQ::INVERSION_OFF) 
      fSUSetSequenceString(pMrProt,pSeqLim,pSeqExpo,"epse");
    else
      fSUSetSequenceString(pMrProt,pSeqLim,pSeqExpo,"epir");        
     
  #endif
     

  #ifdef EP2D_FID
    //fSUSetSequenceString(pMrProt,pSeqLim,pSeqExpo,"epfid");
	pSeqExpo->setSequenceString( "B.Poser@MBIC,NL" );
	pSeqExpo->setSeqShortString( "SMS-EPI" );

    #ifdef PERF
      #ifndef VXWORKS
        pSeqExpo->setApplicationCard     (SEQ::APPLICATION_CARD_PERF);
        pSeqExpo->setApplicationCardName (SEQ::APPLICATION_CARD_NAME_PERF);
      #endif

      //
      //  Due to the new concept of EVA-protocols the perfusion-data was removed from the MrProt.
      //  Now we can not choose between IceProgramOnline2D and IceProgram2D any more.
      //  So we are forced to use the IceProgram2D in any case and have to live with the
      //  following restriction for the moment (or forever?, CHARM 305697):
      //
      //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //  SEQUENCE ep2d_fid CAN ONLY BE EXECUTED FOR IMAGING WITHOUT PERF-POSTPROCESSING WITH
      //  OFFLINE RECONSTRUCTION.
      //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //
      //  This is due to the fact that IceProgram2D does not support online-reconstruction,
      //  but is required to do perfusion postprocessing.
      //
      pSeqExpo->setOnlineFFT         (SEQ::ONLINE_FFT_NONE); 
      pSeqExpo->setICEProgramFilename("%SiemensIceProgs%\\IceProgram2D");

      pSeqExpo->setICEProgramParam   (ICE_PROGRAM_PARA_SHOW_OFFLINE, SEQ::SO_SHOW_YES);
      pSeqExpo->setICEProgramParam   (ICE_PROGRAM_PERF_THRESHOLD  ,    50);
      pSeqExpo->setICEProgramParam   (ICE_PROGRAM_PERF_BASE_START ,     5);
      pSeqExpo->setICEProgramParam   (ICE_PROGRAM_PERF_BASE_END   ,     9);
      pSeqExpo->setICEProgramParam   (ICE_PROGRAM_PERF_TTP_SCALING,100000);
      pSeqExpo->setICEProgramParam   (ICE_PROGRAM_PERF_TTP_START  ,  1000);
    #endif

    #if (defined BOLD) || (defined PASL)
      #ifndef VXWORKS
        pSeqExpo->setApplicationCard     (SEQ::APPLICATION_CARD_FMRI     );
        pSeqExpo->setApplicationCardName (SEQ::APPLICATION_CARD_NAME_BOLD);
	  #endif

      //
      //  Due to the new concept of EVA-protocols the fMRI-data was removed from the MrProt.
      //  Now we can not choose between IceProgramOnline2D and IceProgram2D any more.
      //  So we are forced to use the IceProgramOnline2D in any case and have to live with the
      //  following restriction for the moment (or forever?, see CHARM 305697):
      //
      //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //  SEQUENCE ep2d_bold CAN ONLY BE EXECUTED FOR IMAGING WITHOUT fMRI-POSTPROCESSING WHEN
      //  ONLINE RECONSTRUCTION IS POSSIBLE.
      //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //
      //  This is due to the fact that IceProgramOnline2D does not support OFFLINE-reconstruction,
      //  but is required to do fMRI postprocessing.
      //
      pSeqExpo->setOnlineFFT          (SEQ::ONLINE_FFT_PHASE);
      // JR 4 lines commented out 
	  // #ifdef PACE3D
          pSeqExpo->setICEProgramFilename ("%SiemensIceProgs%\\IceProgram2D");
      //  #else
      //    pSeqExpo->setICEProgramFilename ("%SiemensIceProgs%\\IceProgramOnline2D");
      //#endif

	  //benpos WIP
	  //call modified IceProgram in any case (to support features like Mosaic 'save-uncombines images', mag/phs reconstruction and superlarge matrix size)
	  //benpos caipi -->
	  if ( l_LongEntry_SmsFactor > 1  && b_Checkbox_SmsOnlineRecon )
	  {
	 	pSeqExpo->setICEProgramFilename ("%CustomerIceProgs%\\BPIceProgramMultiEchoMultiSlice");
	  }
	  else 
	  {
	 	pSeqExpo->setICEProgramFilename ("%CustomerIceProgs%\\BPIceProgramMultiEcho_WIP");
	  }


      #ifndef VXWORKS
        if (pSeqExpo->getOnlineFFT()!=SEQ::ONLINE_FFT_PHASE)
        {
            // mySeqLoop.setIceProgram says offline!
            // we can't help, if IceProgramOnline2D supports only online,
            // ergo:
            //
            if (!pSeqLim->isContextPrepForBinarySearch()) textTR("mySeqLoop.setIceProgram says offline");
        
            if (mySeqLoop.getCTRL_SeqLoop_SwitchToOffline())
            {
                textTR("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                textTR("mySeqLoop.setIceProgram says offline, because it is forced to do so");
                textTR("check registry-key: SOFTWARE\\Siemens\\Numaris4\\Config\\Modality\\Sequence\\CTRL_SeqLoop_SwitchToOffline");
                textTR("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            }

            lStatus=SEQU_ERROR;
            goto FINISHED;
        }
      #endif              
      
    #endif
  #endif


  #ifdef EP2D_SE
    if(pMrProt->preparationPulses().inversion() == SEQ::INVERSION_OFF) fSUSetSequenceString(pMrProt,pSeqLim,pSeqExpo,"epse");
    else                                                               fSUSetSequenceString(pMrProt,pSeqLim,pSeqExpo,"epir");

    if (REOInfo.getPhasePartialFourierFactor()<0.70)
    {
        pSeqExpo->setPCAlgorithm(SEQ::PC_ALGORITHM_MARGOSIAN);
    }

    #ifndef VXWORKS
      pSeqExpo->setApplicationCard(SEQ::APPLICATION_CARD_INLINE);
    #endif
  #endif
  
  // ---------------------------------------------------------------------------
  // specify submatrix for PC algorithm
  // ---------------------------------------------------------------------------
  if (pSeqExpo->getPCAlgorithm() != SEQ::PC_ALGORITHM_NONE)
  {
    pSeqExpo->setNoOfPhaseCorrLines   (  16 * pMrProt->kSpace().baseResolution() / 256 );
    pSeqExpo->setLinSlopeLength       (  16 * pMrProt->kSpace().baseResolution() / 256 );
    pSeqExpo->setNoOfPhaseCorrColumns ( 128 * pMrProt->kSpace().baseResolution() / 256 );
    pSeqExpo->setColSlopeLength       ( 128 * pMrProt->kSpace().baseResolution() / 256 );
  }

  // ---------------------------------------------------------------------------  
  // check hard limits for maximum number of receiver channels
  // ---------------------------------------------------------------------------
  if ( configInfo.receiverChannelLimitsSet() )
  {
	// find channel limit for current resolution
	long lMaxChannels = configInfo.getMaxReceiverChannels( pMrProt->kSpace().baseResolution() );
  	      	
  	// too many channels set ? -> cannot do it
  	pSeqExpo->setMaxReceiverChannels( lMaxChannels );
  	if ( pMrProt->coilInfo().Meas().getNumOfUsedRxChan() > lMaxChannels )
  	{
      lStatus=SEQU_ERROR;

	  if ( !pSeqLim->isContextPrepForBinarySearch() )
	  {
		  textTR("maximum number of receiver channels exceeded");                   
	  }

      goto FINISHED;
    }        
  }

  if (lDebug_SEQ_fSEQPrep==39) longTR("pSeqExpo->getMaxReceiverChannels()",pSeqExpo->getMaxReceiverChannels());

#ifdef SUPPORT_iPAT_a_ep2d

  //---------------------------------------------------------------------------
  // call fPATPrepPost, which checks some PAT related restrictions
  //---------------------------------------------------------------------------
  lStatus = fPATPrepPost(pMrProt,pSeqLim,pSeqExpo,&REOInfo);

    #ifdef BOLD 
        // Forbid channel reduction for GRAPPA
        if ( pMrProt->PAT().PATMode() == SEQ::PAT_MODE_GRAPPA )
        {
            // Read ICE control flag
            long lIceControl = pSeqExpo->getICEProgramParam(ICE_PROGRAM_PARA_CTRL_MASK);

            // Set the "no channel reduction" bit
            lIceControl = lIceControl | ICE_PROGRAM_MSK_iPAT_NO_CHANNEL_REDUCTION;
      
            // Write back
            pSeqExpo->setICEProgramParam(ICE_PROGRAM_PARA_CTRL_MASK, lIceControl);
        }
    #endif  // endif BOLD
  
    if (! pSeqLim->isContextPrepForBinarySearch() )  
    {
        CheckStatusPB (lStatus, "fPATPrepPost");
    } 
    else  
    {
        CheckStatusB  (lStatus);
    }
#endif

  // ---------------------------------------------------------------------------
  // Check the required TE, TI and TR values against the current ones, adapt
  // those values, if pSeqLim->isContextPrepForMrProtUpdate() is true.
  // ---------------------------------------------------------------------------
  if (! fEPIStdUICheckNeededTETITR(pMrProt,pSeqLim,lNeededTE,lNeededTI,lNeededTR))
  {
      lStatus = SEQU_ERROR; // no specific error-message available,
                            // not nice but SEQU__NEGATIVE_TE_FILL would also be missleading
  }

  
  
  // ---------------------------------------------------------------------------
  //benpos WIP
  //update parameters in the special card
	 //make sure the delay time variables and their visibility make sense
    
	// register desired combination mode to MrProt so ICE knows it:
	pMrProt->wipMemBlock().alFree[20]  = l_SelectionBox_EchoCombineMode;


	// we cannot use online combination for MAG_AND_PHA recon
	if (pMrProt->reconstructionMode() ==  SEQ::RECONMODE_MAGN_PHASE )
	{
		l_SelectionBox_EchoCombineMode = 0;
	}

	// we cannot use online combination for uncombined (single coil) images 
    if ( pMrProt->uncombImages() )  
	{
		l_SelectionBox_EchoCombineMode = 0;
	}


	//make sure the visibility of variables makes sense
    SHOW_PARAM(pSeqLim, &l_LongEntry_delayBetweenFirstAndSecondEcho, true);
    SHOW_PARAM(pSeqLim, &l_LongEntry_delayBetweenRemainingEchoes, true); 
	SHOW_PARAM(pSeqLim, &l_SelectionBox_EchoCombineMode, true); 
	 

	// set forwarding of individual echoes:  1) user-selectable in case of a ME-EPI combination mode, 2) always true it only one echo (else case),
	if (pMrProt->contrasts() > 1 ) 
	{
		SHOW_PARAM(pSeqLim, &b_Checkbox_SaveSeparateEchoes, true);
		// force to true if no combinatoin mode is selected
		if (l_SelectionBox_EchoCombineMode == 0 ) 
		{
			b_Checkbox_SaveSeparateEchoes = true;
			pMrProt->wipMemBlock().alFree[22]  = 1;
		}
		else // user choice if combination is selected
		{
			if ( b_Checkbox_SaveSeparateEchoes == true  ) pMrProt->wipMemBlock().alFree[22]  = 1;
			if ( b_Checkbox_SaveSeparateEchoes == false ) pMrProt->wipMemBlock().alFree[22]  = 0;
		}
	}
	else //always forward 'uncombined echoes' data in case of only one echo
	{		
		SHOW_PARAM(pSeqLim, &b_Checkbox_SaveSeparateEchoes, false);
		b_Checkbox_SaveSeparateEchoes = true;
		pMrProt->wipMemBlock().alFree[22]  = 1;
	}

	
	// for fixed T2* weighting show the coresponding parameter field
	if (l_SelectionBox_EchoCombineMode == 2 )
	{
		SHOW_PARAM(pSeqLim, &l_LongEntry_T2starweight, true); 
        pMrProt->wipMemBlock().alFree[21]  = l_LongEntry_T2starweight;
	}
	else
	{
		HIDE_PARAM(pSeqLim, &l_LongEntry_T2starweight);
		pMrProt->wipMemBlock().alFree[21]  = 0;
	}

	if (pMrProt->contrasts() < 2 ) 
	{
		l_LongEntry_delayBetweenFirstAndSecondEcho = 0;
		 HIDE_PARAM(pSeqLim, &l_LongEntry_delayBetweenFirstAndSecondEcho);

		 HIDE_PARAM(pSeqLim, &l_SelectionBox_EchoCombineMode);
		 pMrProt->wipMemBlock().alFree[20]  = 0;
		 pMrProt->wipMemBlock().alFree[21]  = 0;
	}
	
	if (pMrProt->contrasts() < 3 ) 
	{
		l_LongEntry_delayBetweenRemainingEchoes    = 0;
		 HIDE_PARAM(pSeqLim, &l_LongEntry_delayBetweenRemainingEchoes);
	}


	SHOW_PARAM(pSeqLim, &l_LongEntry_FixedZshim, true); 
	SHOW_PARAM(pSeqLim, &l_LongEntry_AlternZshim, true); 
	if (pMrProt->contrasts() > 1   ) 
	{
		SHOW_PARAM(pSeqLim, &l_LongEntry_FixedZshim, false); 
		SHOW_PARAM(pSeqLim, &l_LongEntry_AlternZshim, false); 
		l_LongEntry_FixedZshim  = 0;
		l_LongEntry_AlternZshim = 0;
	}









	SHOW_PARAM(pSeqLim, &l_SelectionBox_RefScanMode, false);
	SHOW_PARAM(pSeqLim, &l_LongEntry_FleetDummies, false); 
	SHOW_PARAM(pSeqLim, &l_LongEntry_FleetFlipAngle, false); 
//	SHOW_PARAM(pSeqLim, &d_DoubleEntry_SSBPatRefscanTE, false); 
//	SHOW_PARAM(pSeqLim, &d_DoubleEntry_SSBPatRefscanBW, false); 
//	SHOW_PARAM(pSeqLim, &d_DoubleEntry_SSBPatRefscanBaseRes, false); 
	if (pMrProt->PAT().PATMode() != SEQ::PAT_MODE_NONE )
	{
		SHOW_PARAM(pSeqLim, &l_SelectionBox_RefScanMode, true);
		
		if ( l_SelectionBox_RefScanMode == 0 || l_SelectionBox_RefScanMode == 1 )	// normal segemnted or single shot
		{
			//do nothing secial
		}
		else if (l_SelectionBox_RefScanMode == 2 )									//FLEET
		{ 
			SHOW_PARAM(pSeqLim, &l_LongEntry_FleetDummies, true); 
			SHOW_PARAM(pSeqLim, &l_LongEntry_FleetFlipAngle, true); 
		}		
		else if (l_SelectionBox_RefScanMode == 3 )									// FLASH
		{
	//		pMrProt->wipMemBlock().alFree[8]  = 3;
		//	SHOW_PARAM(pSeqLim, &d_DoubleEntry_SSBPatRefscanTE, true); 
		//	SHOW_PARAM(pSeqLim, &d_DoubleEntry_SSBPatRefscanBW, true); 
		//	SHOW_PARAM(pSeqLim, &d_DoubleEntry_SSBPatRefscanBaseRes, true)
		}
		else { }
	}



	
	SHOW_PARAM(pSeqLim, &l_LongEntry_FatSatFlipAngle, false);
	if (pMrProt->preparationPulses().fatSuppression() == SEQ::FAT_SATURATION ) //the CSatFat I think only does the normal fat sat so the flip angle is only for that 
	{
		SHOW_PARAM(pSeqLim, &l_LongEntry_FatSatFlipAngle, true);	
	}


	// benpos SMS-CAIPI -->
	//sort out parameter visibility for sms 

	SHOW_PARAM(pSeqLim, &l_LongEntry_SmsFactor, false);
	SHOW_PARAM(pSeqLim, &l_LongEntry_SmsShiftFactor, false);
	SHOW_PARAM(pSeqLim, &b_Checkbox_SmsPhaseOptimise, false);
	SHOW_PARAM(pSeqLim, &b_Checkbox_SmsOnlineRecon, false);

	if (pMrProt->sliceSeries().size() > 1 ) 
	{
		SHOW_PARAM(pSeqLim, &l_LongEntry_SmsFactor, true); 
		//TODO: tooltip or read-only parameter for slice separation	
		double SMS_slicegap = pMrProt->sliceSeries().size() / l_LongEntry_SmsFactor * pMrProt->sliceGroupList()[0].distance();
		if ( l_LongEntry_SmsFactor > 1)
		{	
			SHOW_PARAM(pSeqLim, &l_LongEntry_SmsShiftFactor, true);
			SHOW_PARAM(pSeqLim, &b_Checkbox_SmsPhaseOptimise, true);
			SHOW_PARAM(pSeqLim, &b_Checkbox_SmsOnlineRecon, true);
		}
		
		
		//enforce a multiplexing factor that is doable	
		if ( fabs( pMrProt->sliceSeries().size()  % l_LongEntry_SmsFactor ) > 0  )
		{
			do {
				if (l_LongEntry_SmsFactor > 1 ) 	l_LongEntry_SmsFactor=l_LongEntry_SmsFactor  -1 ;
			} 	while (fabs( pMrProt->sliceSeries().size()  % l_LongEntry_SmsFactor ) > 0 );
		}
		
		// TODO: also enforce an odd number of excited slices
		if  ( (pMrProt->sliceSeries().size() / l_LongEntry_SmsFactor)%2 ==0 ) cout << " excited number of slices is EVEN - not nice!!!! " << endl;
		
		
		
		//finally, write SmsFactor and SmsShiftFactor to WipMemBlock
		pMrProt->wipMemBlock().adFree[1]  = double (l_LongEntry_SmsFactor);
		pMrProt->wipMemBlock().adFree[2]  = double (l_LongEntry_SmsShiftFactor);
	
	}

 
	// <-- benpos SMS-CAIPI


	//set numbers of inversions per TR depending on VASO mode
	if (pMrProt->Asl().Mode()==SEQ::ASL_CUSTOM1 ) 
	{
		l_LongEntry_VolumesPerInversion = 1;
	}
	else if (pMrProt->Asl().Mode()==SEQ::ASL_CUSTOM2 )  
	{
		l_LongEntry_VolumesPerInversion = 2;
	}
	else {}

	
	


	pSeqExpo->setAdditionalScaleFactor(d_DoubleEntry_FFTscale);

//////// this here would be for forcing a phase resolution that is larger than baseresolution
/*	
pMrProt->kSpace().phaseEncodingLines(  2* pMrProt->kSpace().baseResolution()    );
cout << " pMrProt->kSpace().baseResolution() " << pMrProt->kSpace().baseResolution() << " pelines() " << pMrProt->kSpace().phaseEncodingLines() << endl;
//pMrProt->sliceSeries()[0].phaseFOV();// 2* pMrProt->sliceSeries().at(0).readoutFOV() ) ;
cout << " pMrProt->sliceSeries()[0].phaseFOV( " << pMrProt->sliceSeries()[0].phaseFOV() << endl;
/*
//pMrProt->sliceSeries().asSlice[0].dPhaseFOV(2*pMrProt->sliceSeries().asSlice[0].dReadoutFOV() ) ;
pMrProt->sliceSeries().slice[0].dPhaseFOV(2*pMrProt->sliceSeries().slice[0].dReadoutFOV() ) ;
*/  ////////////////////////////////////////

	// ggg Beginning of actual spiral code
	gSpiralPrep(pMrProt,pSeqLim,pSeqExpo);
	if(SpiralModeOn()) {
		pMrProt->kSpace().trajectory(SEQ::TRAJECTORY_SPIRAL);
	} else {
		pMrProt->kSpace().trajectory(SEQ::TRAJECTORY_CARTESIAN);
	}

	// run the protocol update
	 UPDATE_PROTOCOL(pMrProt, pSeqLim);

  // ---------------------------------------------------------------------------
  
  // benpos write the MrProt to file
  #ifndef VXWORKS	
	  //sprintf(protfilename, "%s\\RF_pulses\\MrProt_ME-EPI.txt", env_custom_seq );	
	  sprintf(protfilename, "%s\\__MrProt_ME-EPI.txt", env_custom_seq );	
	  pMrProt->fwrite(protfilename); 
  #endif
  
  // ---------------------------------------------------------------------------
  // finished
  // ---------------------------------------------------------------------------
  mTRend;
  return(lStatus);

  // ---------------------------------------------------------------------------
  // in the error case: 
  // ---------------------------------------------------------------------------
FINISHED:
  if (NLS_SEVERITY(lStatus) == NLS_SUCCESS) lStatus = SEQU_ERROR;   // make sure that an ERROR code is returned
  fEPIStdResetSolveHandlerControlTETITR();                          // even if TE,TI or TR is changed, we can't help

  if (mIsTRend) entryTR("aborted!");
  return(lStatus);
}

/*[ Function ****************************************************************\
*
* Name        : fSEQCheck
*               
* Description : Checks the real-time sequence for gradient overflows.
*               
* Return      : An NLS status code.
*
\****************************************************************************/

/*] END: */

NLS_STATUS fSEQCheck
(
  MrProt        *pMrProt,           // IMP: user choice parameters 
  SeqLim        *pSeqLim,           // IMP: limits from fSEQInit() 
  SeqExpo       *pSeqExpo,          // IMP: exports from fSEQPrep()
  SEQCheckMode  *  //pSEQCheckMode      // IMP: check mode 
)
{
  #undef  DEBUG_ORIGIN
  #define DEBUG_ORIGIN 0x00000400
  
  static const char *ptModule = {"fSEQCheck"};
  mTRrun;

  NLS_STATUS   lStatus  = SEQU__NORMAL;
  BOOL         bSuccess = FALSE;


#ifdef EP2D_DIFF
#ifndef VXWORKS
  /* In the past, the test case lGrMomentNotRephasedErr of the SeqUT used to be 
     disabled here exclusively for the AC44, since it could not handle the
     asymmetric Maxwell compensation that has been introduced with CHARM 324742. 
     With CHARM 361864, a generalized Maxwell compensation has been introduced
     that is active for all system types. Instead of disabling selected test
     cases (which - depending on the protocol - yields a wild mixture of
     passed and failed test), the Maxwell correction itself is now disabled in 
     the corresponding diffusion module, but only (!) if the unit test is active.
  */

  if (theConfig.getIntDefault("INIAccess", "Dump", 0)) {
    char tDummy[256];
    theConfig.fDumpAll ( "c:\\temp\\foo.txt" /*stdout*/, tDummy) ;
    printf ("Error=%s", tDummy) ;
  }
#endif
#endif


  // ggg to remove:
  if(SpiralModeOn()) {
	return(lStatus);
  }

  //---------------------------------------------------------------------------
  // Set the looping parameters
  //---------------------------------------------------------------------------
  mySeqLoop.setlinesToCheck        (1);
  mySeqLoop.setlineNoToCheck       (0,0);
  
  mySeqLoop.setpartitionsToCheck   (1);
  mySeqLoop.setparNoToCheck        (0,0); 
  
  //---------------------------------------------------------------------------
  // Execute the check loops
  //---------------------------------------------------------------------------
  bSuccess = mySeqLoop.check
             ( 
                fSEQRunKernel,
                pMrProt,
                pSeqLim,
                pSeqExpo,
                asSLC,
                EPIKernel.getReadOutAddress()
             );   
  if(! bSuccess )mSBBErrGotoFinish( mySeqLoop," mySeqLoop.check")

  //---------------------------------------------------------------------------
  // ready.
  //---------------------------------------------------------------------------
FINISHED:
  mTRend;
  return(lStatus);
}

/*[ Function ****************************************************************\
*
* Name        : fSEQRun
*               
* Description : Executes the real-time sequence.
*               
* Return      : An NLS status code.
*
\****************************************************************************/

/*] END: */

NLS_STATUS fSEQRun
(
  MrProt     *pMrProt,    /* IMP: user choice parameters  */
  SeqLim     *pSeqLim,    /* IMP: limits from fSEQInit()  */
  SeqExpo    *pSeqExpo    /* IMP: exports from fSEQPrep() */
)
{
  #undef  DEBUG_ORIGIN
  #define DEBUG_ORIGIN 0x00002000
  
  static const char *ptModule = {"fSEQRun"};
  mTRrun;

  NLS_STATUS lStatus          = SEQU__NORMAL;
  BOOL       bSuccess         = FALSE;

  // ggg
  gSpiral_Run(pMrProt,pSeqLim,pSeqExpo);
  
  // ---------------------------------------------------------------------------
  // provide unit test with special information:
  // ---------------------------------------------------------------------------
#ifndef VXWORKS
#ifdef EP2D_DIFF
  if ( mIsUnittestActive() ) 
  {
      if (pMrProt->diffusion().mode()==SEQ::DIFFMODE_TENSOR || pMrProt->diffusion().mode()==SEQ::DIFFMODE_FREE) 
      {
          SeqUT.setDiffICEDimensionMode(1,1,0);
          SeqUT.setSizeOfDimFree (0);
          SeqUT.setSizeOfDimSet  (0);
          
          // The sequences uses the Free loop to step through b-values and diffusion 
          // directions. However, in MDDW-mode, this information is exported by using
          // the Crep information of the Mdh. For Crep > 1, SeqUT expects 
          // LastScanInMeas- and LastScanInConcat-flags, which are of course not
          // set by the Free loop. Thus, these test cases have to be disabled. Also,
          // Crep dimension information does not match. 
          //
          // For the future, one should consider a more consistent communication
          // mechanism between sequence and Ice. Optionally, libSeqUT could be
          // adapted to handle this special case appropriately.
          SeqUT.DisableTestCase (lAddDimIndexOutOfRangeErr,  RTEB_ORIGIN_fSEQRunFinish,
              "SeqLoop uses mixes Crep with the free loop counter, thus additional dimension indices do not match");
          SeqUT.DisableTestCase (lMissingLastScanInMeasFlag, RTEB_ORIGIN_fSEQRunFinish,
              "SeqLoop uses mixes Crep with the free loop counter, thus no LastScanInMeas-flags are set");
          SeqUT.DisableTestCase (lNoOfLastScanInConcatErr,   RTEB_ORIGIN_fSEQRunFinish,
              "SeqLoop uses mixes Crep with the free loop counter, thus no LastScanInConcat-flags are set");
          SeqUT.DisableTestCase (lNoAddDimCheckedMomentErr,  RTEB_ORIGIN_fSEQRunKernel,
              "SeqLoop uses mixes Crep with the free loop counter, thus additional dimensions cannot be checked");
      } 
      else 
      {
          SeqUT.setDiffICEDimensionMode(Diff->getNoOfDirections(),
              pMrProt->diffusion().diffWeightings(),
              pMrProt->diffusion().bValue()[0] );
      }
  }  
#endif
#endif
  
  // ---------------------------------------------------------------------------
  // Disable unit test TR calculation when using multiple 
  // concatenations with long TR triggering mode
  // ---------------------------------------------------------------------------
#ifndef VXWORKS

	if ( mySeqLoop.isLongTRTrigMode() && (pMrProt->concatenations() > 1) )
	{
		SeqUT.DisableTestCase(lTRClockErr, RTEB_ClockCheck, "The unit test TR calculations are not valid for multiple concatenations in long TR triggering mode");
	}

#endif

  // ---------------------------------------------------------------------------
  // Enable ADC Packaging
  // ---------------------------------------------------------------------------
  if( pMrProt->NavigatorParam().RespComp() == SEQ::RESP_COMP_OFF || pMrProt->NavigatorParam().RespComp() == SEQ::RESP_COMP_BREATH_HOLD )
  {
      // pack 16 ADCs together for data transfer from PCI-card to-MRIR RAM
 //     fRTSetReadoutPackaging(16); //BENPOS NEED TO DISABLE FOR PTX VOP supervisions
  }

  //---------------------------------------------------------------------------
  // initialize alPrepScanCounter
  //---------------------------------------------------------------------------
  #ifdef SUPPORT_iPAT_a_ep2d
  {
    for (long lK=0; lK<K_NO_SLI_MAX; lK++) {alPrepScanCounter[lK]=0;}
  }
  #endif

  //---------------------------------------------------------------------------
  // initialize alPaslPrepScanCounter, alPaslScanCounter
  //---------------------------------------------------------------------------
  #ifdef PASL
  {
    for (long lK=0; lK <K_NO_SLI_MAX; lK++ ) {alPaslPrepScanCounter[lK]=0;}
    for (long lK2=0;lK2<K_NO_SLI_MAX; lK2++) {alPaslScanCounter[lK2]   =0;}
  }
  #endif

  // ---------------------------------------------------------------------------
  // initialize osc-bit control-flags
  // ---------------------------------------------------------------------------
  #ifndef EP2D_DIFF
  {
    for (long lI=0; lI<lMaxOscBitSentFlags; lI++)
    {
        abOscBitSentForMeas[lI]=false;
    }
  }
  #endif

  //---------------------------------------------------------------------------
  // Initialization of the unit test function
  //---------------------------------------------------------------------------
  stringTR("running", pSeqLim->getLinkedSeqFilename() );
  mSEQTest(pMrProt,pSeqLim,pSeqExpo,RTEB_ORIGIN_fSEQRunStart,0,0,0,0,0); /*! EGA-All !*/

  
   
    
  //---------------------------------------------------------------------------
  // benpos WIP: start logging the physio data, after creating the files with timestamp
  //---------------------------------------------------------------------------
    #ifdef VXWORKS
    
    if (b_Checkbox_LogPhysio)
	{
		time_t rawtime;
		struct tm * timeinfo;
		char timestamp[512];
		char pulsfile[512];
		char respfile[512];
		char ecgfile[512];
		char extfile[512];

		time ( &rawtime );
		timeinfo = gmtime ( &rawtime );
		sprintf(timestamp, "%d%02d%02d_%02d%02d%02d", timeinfo->tm_year+1900, 
													   timeinfo->tm_mon+1,
													   timeinfo->tm_mday,
													   timeinfo->tm_hour,
													   timeinfo->tm_min,
													   timeinfo->tm_sec);
		sprintf(pulsfile, "Physlog\\Pulslog_%s", timestamp);
		sprintf(respfile, "Physlog\\Resplog_%s", timestamp);
		sprintf(ecgfile, "Physlog\\ECGlog_%s", timestamp);
		sprintf(extfile, "Physlog\\EXTlog_%s", timestamp);
		MyPhysioPuls.startLoggingSignal(SEQ::SIGNAL_PULSE, pulsfile);
  		MyPhysioResp.startLoggingSignal(SEQ::SIGNAL_RESPIRATION, respfile);
		MyPhysioECG.startLoggingSignal( SEQ::SIGNAL_EKG, ecgfile);
		MyPhysioEXT.startLoggingSignal( SEQ::SIGNAL_EXT, extfile);
    }  
   #endif

  
  
  gSpiral_setSomeParamsBeforeKernel(1,1);
  gSpiral_setRelevantForMeasTime();

  
  //---------------------------------------------------------------------------
  // Execute the measurement loops
  //---------------------------------------------------------------------------
  bSuccess = mySeqLoop.run_new
             (
                (NLS_STATUS (*)(TYPESFor_fSEQRunKernel)) fSEQRunKernel,
                pMrProt,
                pSeqLim,
                pSeqExpo,
                asSLC,
                EPIKernel.getReadOutAddress()
             );
  if(! bSuccess )mSBBErrGotoFinish( mySeqLoop," mySeqLoop.run")


  //---------------------------------------------------------------------------
  // benpos WIP: stop logging the physio data
  //---------------------------------------------------------------------------

    #ifdef VXWORKS
	if (b_Checkbox_LogPhysio)
	{
		MyPhysioPuls.stopLoggingSignal(SEQ::SIGNAL_PULSE);
		MyPhysioResp.stopLoggingSignal(SEQ::SIGNAL_RESPIRATION);
		MyPhysioECG.stopLoggingSignal(SEQ::SIGNAL_EKG);
		MyPhysioEXT.stopLoggingSignal(SEQ::SIGNAL_EXT);
	} 
    #endif


  //---------------------------------------------------------------------------
  // ready.
  //---------------------------------------------------------------------------
FINISHED:
  mSEQTest(pMrProt,pSeqLim,pSeqExpo,RTEB_ORIGIN_fSEQRunFinish,0,0,0,0,0); /*! EGA-All !*/
  stringTR("finished", pSeqLim->getLinkedSeqFilename() );
  mTRend;
  return(lStatus);
}

/*[ Function ****************************************************************\
*
* Name        : fSEQRunKernel
*               
* Description : Executes the basic timing of the real-time sequence.
*               This function is called by the function (libSBB)fSEQRunStdTseIR.
*               
* Return      : An NLS status code.
*
\****************************************************************************/

/*] END: */
static NLS_STATUS fSEQRunKernel
(
  MrProt  *pMrProt,
  SeqLim  *pSeqLim,
  SeqExpo *pSeqExpo,
  long     lKernelMode,
  long     lSlice,
  //ggg
  long     lPartition,
  long     lLine
  //long     ,             //lPartition
  //long     //lLine
)
{
  #undef  DEBUG_ORIGIN
  #define DEBUG_ORIGIN 0x00004000
  
  static const char *ptModule = {"fSEQRunKernel"};
  mTRrun;

  

  // benpos SBBPatRefScan - first thing: run SBBPatRefScan here
	if ( l_SelectionBox_RefScanMode == 3 && !b_sbbpatrefscans_done )
	{
		fRTSetReadoutEnable(1);
		// loop over all slices and partitions
		for (long lSliceIndex = 0; lSliceIndex < pMrProt->sliceSeries().size(); lSliceIndex++  ) 
		{
			SBBPATRefScan.setSlice(lSliceIndex);
  			//for( long lPartitionsCounter = 0; lPartitionsCounter < SBBPATRefScan.getlPartitionsToMeasure(); lPartitionsCounter++)
			//{
			//	SBBPATRefScan.setPartition(lPartitionsCounter);
				if( !SBBPATRefScan.run(pMrProt, pSeqLim, pSeqExpo, asSLC) )
				{
					if ( ! pSeqLim->isContextPrepForBinarySearch() )
					{
					   TRACE_PUT1(TC_ALWAYS, TF_SEQ,"%s: Error encountered in call of SBBPATRefScan.run()", ptModule );
					}                  
					return  SBBPATRefScan.getNLSStatus() ;
				}
			//}        
		}
		if ( lKernelMode != KERNEL_CHECK  ) 
		{
			b_sbbpatrefscans_done = true;
		}
		
		// wait some time before the first volue (signal is low due to preceeding trufi-readout
		long lWait = 10000;
		fRTEBInit(&asSLC[lSlice].m_sROT_MATRIX);
		fRTEI( lWait , 0, 0, 0, 0, 0, 0, 0); 
		fRTEBFinish();

		fRTSetReadoutEnable(0);
		
	}
// end benpos

  //---------------------------------------------------------------------------
  // to send phase-FT flags correctly:
  //---------------------------------------------------------------------------
  REOInfo.setIsLastAcquisition(EPIKernel.getReadOutAddress()->Mdh.getCacq() == pMrProt->averages()-1);


#ifdef SUPPORT_PACE
  if( pMrProt->NavigatorParam().RespComp() != SEQ::RESP_COMP_OFF )
  {
      EPIKernel.setADCRelevant(mySeqLoop.isFirstADCRelevant(),mySeqLoop.isLastADCRelevant());
  }
#endif
  
  #ifdef EP2D_DIFF  
  //---------------------------------------------------------------------------
  // Set diffusion loop counter:
  //
  // Since seqloop does not perform a loop over the free counter during check 
  // phase we have to set one b-value to calc the stimulation prediction 
  // correctly. If there is more than one b value, we take the last b-value, 
  // which is the highest, and consider it to be the most stimulating one...
  //--------------------------------------------------------------------------

  static long lCurrentRepetition = -1;     

  // Are we using the old ICE program or the new online program?
  bool bOldIceProgram = ( NULL != strstr (pSeqExpo->getICEProgramFilename(), "ProgramDiffusion2D") ) ;


  if ( (lKernelMode & KERNEL_CHECK) == KERNEL_CHECK )
  {
      lDiffLoopCounter  = maximum (0L, (long)Diff->getTotalScans(pMrProt)-1L) ;
  }
  else
  {
      lDiffLoopCounter  = mySeqLoop.getDiffusionLoopCounter() ;

	  if (! bOldIceProgram &&
          (pMrProt->diffusion().mode()==SEQ::DIFFMODE_TENSOR ||
	       pMrProt->diffusion().mode()==SEQ::DIFFMODE_FREE) ) 
	  { 
		      if (lCurrentRepetition==-1)
			    lCurrentRepetition = EPIKernel.getReadOutAddress()->Mdh.getCrep() ;
		      
			  
			  EPIKernel.getReadOutAddress()->Mdh.setCrep (
                 (unsigned short)(lCurrentRepetition*(long)Diff->getTotalScans(pMrProt) + lDiffLoopCounter)) ;
  
	          if (EPIKernel.getReadOutAddress()->Mdh.getEvalInfoMask() & MDH_LASTSCANINMEAS) 
			  {
                  lCurrentRepetition = -1 ;  // enforce a re-read
                  // TRACE_PUT0 (TC_INFO, TF_SEQ, "*** LastScanInMeas") ;
			  }
	  }
  }

  
  mPrintTrace5 (DEBUG_RUN, DEBUG_LANDMARK, 
                "Scan SeqLoop #%ld, Cset=%d, Cfree=%d, Crep=%d    (Kernelmode=%ld)",
                mySeqLoop.getDiffusionLoopCounter(), 
                EPIKernel.getReadOutAddress()->Mdh.getCset(),
                EPIKernel.getReadOutAddress()->Mdh.getCida(),
				EPIKernel.getReadOutAddress()->Mdh.getCrep(),
                lKernelMode ) ;

  #endif


  // ggg
  gSpiral_setRep(EPIKernel.getReadOutAddress()->Mdh.getCrep());
  
  //---------------------------------------------------------------------------
  // send clock check event for unit test (SBBBinomialPulses does not do it)
  //---------------------------------------------------------------------------
  mSEQTest(pMrProt,pSeqLim,pSeqExpo,RTEB_ClockCheck,27,-1,asSLC[lSlice].getSliceIndex(),0,0); /*! EGA-All !*/

  //---------------------------------------------------------------------------
  // remaining kernel configurations
  //---------------------------------------------------------------------------
  #ifdef PASL       
  // disable seqloop TR fill to shift slices together, add lTRfillPASL after last slice 
  static long lTRfillPASL = 0;  
  lTRfillPASL = mySeqLoop.getlTRFill();     
  #else
  EPIKernel.setTRFill(mySeqLoop.getlTRFill());
  #endif

  #ifdef EP2D_DIFF
  // switch off blips for phase corr scan (only if not PAT)
  // if (pMrProt->PAT().PATMode()==SEQ::PAT_MODE_NONE) 
  // EPIKernel.setExecuteKernelAsPhaseCorrectionScan (lKernelMode==KERNEL_PHASECOR);
  #else
  {
    // CHARM 304696: send osc bit only once per measurement
    //
    if (EPIKernel.getReadOutAddress()->Mdh.getCrep() > lMaxOscBitSentFlags-1)
    {
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s: EPIKernel.getReadOutAddress()->Mdh.getCrep() > lMaxOscBitSentFlags-1",ptModule);
        TRACE_PUT1(TC_INFO, TF_SEQ, "%s: => sending osc-bit uncontrolled",ptModule);
    
        EPIKernel.setDoNotSendOscBit    (false);
        EPIKernel.setDoNotSendExtTrigger(false);
    }
    else
    {
        if (abOscBitSentForMeas[EPIKernel.getReadOutAddress()->Mdh.getCrep()] || lKernelMode == KERNEL_PREPARE )
        {
            EPIKernel.setDoNotSendOscBit    (true);
            EPIKernel.setDoNotSendExtTrigger(true);
        }
        else
		{
			//benpos add this extra IF-ELSE loop to make sure we only send ONE trigger per multi-TI readout (in case that is on)
			if (    EPIKernel.getReadOutAddress()->Mdh.getCrep() % l_LongEntry_VolumesPerInversion == 0 )
			{
				EPIKernel.setDoNotSendOscBit    (false);
				EPIKernel.setDoNotSendExtTrigger(false);
				abOscBitSentForMeas[EPIKernel.getReadOutAddress()->Mdh.getCrep()]=true;   
			} //<-benpos
			else
			{
				
				EPIKernel.setDoNotSendOscBit    (true);
				EPIKernel.setDoNotSendExtTrigger(true);
			}//<-benpos
		}
    }
  }
  #endif
  
  //----------------------------------------------------------------------------
  // adapt Crep counter and Mdh flags when using PASL reference scan
  // NOTE: needs to be done BEFORE gPaceFeedback.SyncAndIncorporateFeedback !
  //----------------------------------------------------------------------------
  // if this is the first slice, execute SBBPasl (which includes FatSat+Spoiler)
  //----------------------------------------------------------------------------
  #ifdef PASL   
    if( bEnableFirstPrepScanAsM0Scan )
    {     
      if( lKernelMode == KERNEL_PREPARE )
      {
        if( alPaslPrepScanCounter[lSlice] == 0 )
        {
          //-------------------------------------------------------------------------------------
          // enable readout for very first M0 scan (this is even before PAT: lPrepScansNoPATRefScans)
          //-------------------------------------------------------------------------------------
          fRTSetReadoutEnable(1);
          mySeqLoop.setRepetitionValueForMdh( 0 );
          EPIKernel.getReadOutAddress()->Mdh.setCrep( (unsigned short)(0) );
          if( lSlice >= (pMrProt->sliceSeries().size()-1) )
          {
            EPIKernel.getReadOutAddress()->Mdh.setLastScanInMeas(true);
            EPIKernel.getReadOutAddress()->Mdh.setLastScanInConcat(true);
            if( mySeqLoop.getNumberOfRelevantADCs() > 0)
              EPIKernel.getReadOutAddress()->setRelevantForMeasTime(true);
          }
        }

        alPaslPrepScanCounter[lSlice]++;
      }
      else if( lKernelMode == KERNEL_IMAGE )
      {
        // increase counter first (start at crep=1): first crep is M0 reference scan
        alPaslScanCounter[lSlice]++;
        
        mySeqLoop.setRepetitionValueForMdh( alPaslScanCounter[lSlice] );
        EPIKernel.getReadOutAddress()->Mdh.setCrep( (unsigned short)(alPaslScanCounter[lSlice]) );
      }
    }
  #endif

  //-------------------------------------------------------------------------------------
  // synchonize PACE feedbacks and sequence here; i.e. for every new measurement starting          
  //-------------------------------------------------------------------------------------
  #ifdef PACE3D
    if ( !gPaceFeedback.SyncAndIncorporateFeedback(pMrProt, pSeqLim, pSeqExpo, asSLC, &EPIKernel, lKernelMode) )
    {
      textTR("PaceFeedback.SyncAndIncorporateFeedback failed.");
      return SEQU_ERROR;
    }
  #endif


/*
  //-------------------------------------------------------------------------------------
  // execute PASL module before first slice and handle alternate label/control state
  //-------------------------------------------------------------------------------------
  #ifdef PASL   
    if( lSlice == 0 ) 
    {
      // set PASL to label/control scan -> first is M0 ASL_REFERENCE, 2nd is "label", 3rd is "control" etc etc
      //    in PostProc: ASL_LABEL scan is subtracted from ASL_CONTROL (3-2, 5-4, 7-6 ...)
      if( lKernelMode == KERNEL_IMAGE )
      {
        if( ((EPIKernel.getReadOutAddress()->Mdh.getCrep() % 2) != 0 ) || pMrProt->repetitions() == 0)
        {
          Pasl.setLabelState( ASLSTATE::ASL_LABEL );
        }
        else
        {
          Pasl.setLabelState( ASLSTATE::ASL_CONTROL );
        }

		
      }
 
      else    // PrepScans
      {
	    		// only a single reference scan
        if( bEnableFirstPrepScanAsM0Scan && alPaslPrepScanCounter[lSlice] <= 1)
        {     
          Pasl.setLabelState( ASLSTATE::ASL_REFERENCE );
        }
        else 
        {
          Pasl.setLabelState( ASLSTATE::ASL_CONTROL );
        }
	  }
      
	  if( !Pasl.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice]) )
      {
        textTR("Pasl.run failed.");
        return (Pasl.getNLSStatus());
      }

      

      #if (defined DEBUG) && (defined PASL_DEBUG_COUT)
        if(EPIKernel.getReadOutAddress()->Mdh.getCrep() == 1)
          cout << Pasl;
      #endif
    }
	else
	{// FATSAT_ALL_SLICES
      // run fatsat for each slice EXCEPT the first slice, which is included in the PASL SBB
      if (pMrProt->preparationPulses().fatSatMode() == SEQ::FAT_SAT_STRONG)
	  {
		  if(! CSatFat.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice])   )
		  {
			  textTR("CSatFat.run failed.");
			  return SEQU_ERROR; 
		  }
		  if(! SpoilGrad.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice]) )
		  {
			  textTR("SpoilGrad.run failed.");
			  return SEQU_ERROR; 
		  }
	  }
	}            
  #endif  //PASL
*/
  
// benpos SMS-CAIPI -->
  // set the default values for the slice-accelerated runs these are modified for different prep scans further down
    if (l_LongEntry_SmsFactor > 1.0)
	{
      EPIKernel.setSmsRunMode(true);

#ifdef EP2D_DIFF
      Diff->setSmsRunMode(true);
#endif 
        mySeqLoop.setInnerSliceNumber(pMrProt->sliceSeries().size()/l_LongEntry_SmsFactor );
		// EPIKernel.setFlagPCforRTFeedback( m_bB0Correction );
//  }
//
//  if (l_LongEntry_SmsFactor > 1.0)
//  {
      if ( ( lKernelMode  == KERNEL_PREPARE ) &&     ( alPrepScanCounter[lSlice] < m_lInitialDummyScans  )  ) 
      {
          EPIKernel.setSmsRunMode(false);
#ifdef EP2D_DIFF
          Diff->m_bDiffusionGradientsEnabled = false; 
          Diff->setSmsRunMode(false);
#endif 
          mySeqLoop.setInnerSliceNumber(pMrProt->sliceSeries().size());
		 // EPIKernel.setFlagPCforRTFeedback( false );
          //the 200 fill time for calibration scans is required to have a 200 us gap between ADC and RF if SPAIR is selected
          EPIKernel.setTRFill(200);
      }
      else if ( ( lKernelMode  == KERNEL_PREPARE ) &&    ( alPrepScanCounter[lSlice]  == m_lInitialDummyScans )  ) 
      {
          // Enable readout
          fRTSetReadoutEnable(1);
          EPIKernel.setSmsRunMode(false);
#ifdef EP2D_DIFF
          Diff->m_bDiffusionGradientsEnabled = false; 
          Diff->setSmsRunMode(false);
#endif 
          mySeqLoop.setInnerSliceNumber(pMrProt->sliceSeries().size());
		  // EPIKernel.setFlagPCforRTFeedback( false );
          EPIKernel.setTRFill(200);
      }
  }
 // <-- benpos SMS-CAIPI

  //-------------------------------------------------------------------------------------
  // Settings for PAT reference scans.
  // For ep2d_diff the b-value is set to zero for PAT ref scans
  //-------------------------------------------------------------------------------------
  
  
  // benpos SMS-CAIPI --> various small modifications here

  #ifdef SUPPORT_iPAT_a_ep2d

    long lFirstPATRefScan = m_lInitialDummyScans;
    long lLastPATRefScan  = lFirstPATRefScan + m_lPATPrepScans - 1;
	EPIKernel.setRunModeFLEET(false);

    if (l_LongEntry_SmsFactor > 1.0)
	{
	  lFirstPATRefScan += m_lSliceAccPrepScans;
	  lLastPATRefScan += m_lSliceAccPrepScans;
	}
	long lFirstPATRefScanExcitation = lFirstPATRefScan * pMrProt->sliceSeries().size() ; 
	long lLastPATRefScanExcitation  = (lLastPATRefScan+1) * pMrProt->sliceSeries().size() - 1 ; 

//cout <<  "kernel mode "  << lKernelMode<< " lPrepScanExcitationCounter " <<lPrepScanExcitationCounter <<  "  lFirstPATRefScanExcitation " << lFirstPATRefScanExcitation << " lLastPATRefScanExcitation " << lLastPATRefScanExcitation << " NPatRefExc (last- first + 1 ) " << lLastPATRefScanExcitation - lFirstPATRefScanExcitation + 1 << endl;

    if (REOInfo.isPATActive() && REOInfo.getPATAccelerationFactorPE()>1)
    {
		
		
		// JB ME-SMS combined
       // if (lKernelMode==KERNEL_PREPARE && alPrepScanCounter[lSlice]>=lPrepScansNoPATRefScans) //benpos SMS-CAIPI 
          if ( ( lKernelMode                 == KERNEL_PREPARE   ) && 
              // ( alPrepScanCounter[lSlice] >= lFirstPATRefScan ) &&
              // ( alPrepScanCounter[lSlice] <= lLastPATRefScan  )
			   ( lPrepScanExcitationCounter >= lFirstPATRefScanExcitation  ) &&
               ( lPrepScanExcitationCounter <= lLastPATRefScanExcitation  )
                 )
        {
            //for regular EPI refscans
			if (l_SelectionBox_RefScanMode < 2 ) 
			{

				fRTSetReadoutEnable(1);

				//long lPATRefScanCounterInSegment = alPrepScanCounter[lSlice]-lPrepScansNoPATRefScans; //benpos SMS-CAIPI 
				  long lPATRefScanCounterInSegment = alPrepScanCounter[lSlice] - lFirstPATRefScan;

				REOInfo.setPATReorderIndexOffsetForRefScans(lPATRefScanCounterInSegment);
				EPIKernel.setCounterInSegmentForEchoShifting(lPATRefScanCounterInSegment);
				EPIKernel.setBlindImagingADCs(REOInfo.getPATBlindADCsBeforeRefScans(), REOInfo.getPATBlindADCsAfterRefScans());
        
//cout << " patref slice = " << lSlice << " seg " <<  lPATRefScanCounterInSegment<<  " alPrepScanCounter[lSlice] " << alPrepScanCounter[lSlice] <<endl;


				  #ifdef PASL
				  //The reference scan must have repetition 0
				  EPIKernel.getReadOutAddress()->Mdh.setCrep ( 0 );
				  #endif

				// select appropriate diffusion SBB for PAT reference scans
				#ifdef EP2D_DIFF
					Diff->m_bDiffusionGradientsEnabled = false ;
				   // Repetitions loop counter of iPAT reference scans is always zero
				  // (setLoopCounters might have set a different value)
				  EPIKernel.getReadOutAddress()->Mdh.setCrep ( 0 );
				#endif

				if (l_LongEntry_SmsFactor > 1.0)
				{
				#ifdef EP2D_DIFF
					Diff->m_bDiffusionGradientsEnabled = false; 
					Diff->setSmsRunMode(false);
				#endif
					
					EPIKernel.setSmsRunMode(false);
					mySeqLoop.setInnerSliceNumber(pMrProt->sliceSeries().size());
					//	EPIKernel.setFlagPCforRTFeedback( false );
					EPIKernel.setTRFill(200);
				}

			}
			else if (l_SelectionBox_RefScanMode == 2 )  //FLEET refscans
			{
				lSlice  = (lPrepScanExcitationCounter-lFirstPATRefScanExcitation) / ( l_LongEntry_FleetDummies + REOInfo.getPATAccelerationFactorPE()   ) ;
				EPIKernel.getReadOutAddress()->Mdh.setCslc(lSlice);
											
				//flag kernel to run the low flipangle pulse.
				EPIKernel.setRunModeFLEET(true);


				long lPATRefScanCounterInSegment = 0;
				
				if (alPrepScanCounter[lSlice] <l_LongEntry_FleetDummies+lFirstPATRefScan )
				{
					fRTSetReadoutEnable(0);
//cout << " dummy slice = " << lSlice << endl;
				}
				else
				{
					fRTSetReadoutEnable(1);
					lPATRefScanCounterInSegment = ( alPrepScanCounter[lSlice] - l_LongEntry_FleetDummies - lFirstPATRefScan ) % REOInfo.getPATAccelerationFactorPE(); 
//cout << " fleet slice = " << lSlice << " seg " <<  lPATRefScanCounterInSegment<<  endl;
				}



				REOInfo.setPATReorderIndexOffsetForRefScans(lPATRefScanCounterInSegment);
				EPIKernel.setCounterInSegmentForEchoShifting(lPATRefScanCounterInSegment);
				EPIKernel.setBlindImagingADCs(REOInfo.getPATBlindADCsBeforeRefScans(), REOInfo.getPATBlindADCsAfterRefScans());
        
				  #ifdef PASL
				  //The reference scan must have repetition 0
				  EPIKernel.getReadOutAddress()->Mdh.setCrep ( 0 );
				  #endif

				// select appropriate diffusion SBB for PAT reference scans
				#ifdef EP2D_DIFF
					Diff->m_bDiffusionGradientsEnabled = false ;
				   // Repetitions loop counter of iPAT reference scans is always zero
				  // (setLoopCounters might have set a different value)
				  EPIKernel.getReadOutAddress()->Mdh.setCrep ( 0 );
				#endif

				if (l_LongEntry_SmsFactor > 1.0)
				{
				#ifdef EP2D_DIFF
					Diff->m_bDiffusionGradientsEnabled = false; 
					Diff->setSmsRunMode(false);
				#endif
					
					EPIKernel.setSmsRunMode(false);
					mySeqLoop.setInnerSliceNumber(pMrProt->sliceSeries().size());
					//	EPIKernel.setFlagPCforRTFeedback( false );
					EPIKernel.setTRFill(200);

				}


								
			}
			else 
			{ }

        }
        else
        {
            // Why we don't use 
            // EPIKernel.setCounterInSegmentForEchoShifting(REOInfo.getPATRefCounterInSegmentWithKSCenter());
            // in the following line is explained where EPIKernel.setUseEchoShifting is used in fSEQPrep.
            //
            REOInfo.setPATReorderIndexOffsetForImagingScans();
            EPIKernel.setCounterInSegmentForEchoShifting(0);
            EPIKernel.setBlindImagingADCs(0,0);

			// select standard diffusion SBB
		    #ifdef EP2D_DIFF
				Diff->m_bDiffusionGradientsEnabled = true ;
            #endif
        }
    }
    else
    {
        REOInfo.setPATReorderIndexOffsetForImagingScans();
        EPIKernel.setCounterInSegmentForEchoShifting(0);
        EPIKernel.setBlindImagingADCs(0,0);

		// select standard diffusion SBB
		#ifdef EP2D_DIFF
			Diff->m_bDiffusionGradientsEnabled = true ;
        #endif
    }

  #endif

  
  //-------------------------------------------------------------------------------------
  // benpos ASL --> insert this here: execute PASL module before first slice and handle alternate label/control state
  //-------------------------------------------------------------------------------------
  #ifdef PASL   


//cout << " m_lInitialDummyScans " << m_lInitialDummyScans << " m_lSliceAccPrepScans " << m_lSliceAccPrepScans << " m_lPATPrepScans " << m_lPATPrepScans << " m_lPostDummyScans " << m_lPostDummyScans << endl;

    if (   lSlice == 0 )
	{
      // set PASL to label/control scan -> first is M0 ASL_REFERENCE, 2nd is "label", 3rd is "control" etc etc
      //    in PostProc: ASL_LABEL scan is subtracted from ASL_CONTROL (3-2, 5-4, 7-6 ...)
      if( lKernelMode == KERNEL_IMAGE )
      {
        
		 //regular case, one volume TR
		if (l_LongEntry_VolumesPerInversion == 1  ) 
		{
			if( ((EPIKernel.getReadOutAddress()->Mdh.getCrep() % 2) != 0 ) || pMrProt->repetitions() == 0)
			{
          //Pasl.setLabelState( ASLSTATE::ASL_REFERENCE );
		  Pasl.setLabelState( ASLSTATE::ASL_CONTROL );// renzo
			}
			else
			{
			  Pasl.setLabelState( ASLSTATE::ASL_CONTROL );
			}
		}
		else // we have multipe Reps per Inversion
		{
			if( ( (EPIKernel.getReadOutAddress()->Mdh.getCrep() /l_LongEntry_VolumesPerInversion)  % 2) != 0  	)
			{
          //Pasl.setLabelState( ASLSTATE::ASL_REFERENCE );
		  Pasl.setLabelState( ASLSTATE::ASL_CONTROL );// renzo
			}
			else
			{
			  Pasl.setLabelState( ASLSTATE::ASL_CONTROL );
			}
		}

		//special case: if this is VASO, we always control
		if ( pMrProt->Asl().Mode() == SEQ::ASL_CUSTOM1 || pMrProt->Asl().Mode() == SEQ::ASL_CUSTOM2 ) 
		{
		  Pasl.setLabelState( ASLSTATE::ASL_CONTROL );
		}
				
      }
 
      else    // PrepScans
      {
	    // only a single reference scan
        if( bEnableFirstPrepScanAsM0Scan && alPaslPrepScanCounter[lSlice] <= 1)
        {     
          //Pasl.setLabelState( ASLSTATE::ASL_REFERENCE );
		  Pasl.setLabelState( ASLSTATE::ASL_CONTROL );// renzo
        }
		
        else 
        {
          Pasl.setLabelState( ASLSTATE::ASL_CONTROL );
        }

	  }

	  
      
	    
	  //benpos - only label for image scans, PostDummyScans, or for InitialDummyScans in case both iPAT abd SliceAcc are off pppppppppp
	  // Gilad reomve ASL
	  if(!TIMode()) {
		  if (  
			     (lKernelMode == KERNEL_IMAGE )
			  ||  (l_LongEntry_SmsFactor == 1 &&  REOInfo.getPATAccelerationFactorPE() == 1 )
			  || ( alPrepScanCounter[lSlice]  > m_lInitialDummyScans + m_lSliceAccPrepScans + m_lPATPrepScans - 1  )
			 )
		  {
	 
			  //benpos multi-TI  - excecute as normal if there is only 1 volume per inversion
			  if (l_LongEntry_VolumesPerInversion == 1 || lKernelMode == KERNEL_PREPARE  ) 
			  {
				  if (  EPIKernel.getReadOutAddress()->Mdh.getCrep() %2 == 0) { // renzo new loop 
					  if( !Pasl.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice]) )
					  {
					  textTR("Pasl.run failed.");
						return (Pasl.getNLSStatus());
					}
				  }
			  }
			  else // l_LongEntry_VolumesPerInversion > 1 --> run only when the time has come
			  {
				  if ( EPIKernel.getReadOutAddress()->Mdh.getCrep()%l_LongEntry_VolumesPerInversion == 0 )
				  {
					  if( !Pasl.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice]) )
					  {
					    textTR("Pasl.run failed.");
					return (Pasl.getNLSStatus());
					  }
				  }
					
			  }
		  }
		} // ASL only if not TIMode()
    }
	// Gilad run FatSat also on first
	// else 
	if( (lSlice != 0) || TIMode() )
	{//  -- not the first slice:" play FatSat if FatSat Mode is STRONG 
      // run fatsat for each slice EXCEPT the first slice, which is included in the PASL SBB
      if (pMrProt->preparationPulses().fatSatMode() == SEQ::FAT_SAT_STRONG)
	  {
		  if(! CSatFat.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice])   )
		  {
			  textTR("CSatFat.run failed.");
			  return SEQU_ERROR; 
		  }
		  if(! SpoilGrad.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice]) )
		  {
			  textTR("SpoilGrad.run failed.");
			  return SEQU_ERROR; 
		  }
	  }
	}    
	
/* / HERE WE TRY TO RUN FATSAT FOR THE PREPSCANS
// problems: 
//  - FatSat is only prepper in case of STRONG
//  - then Energy gets automatically added for each slice (each kernelcall) also for the time series


//benpos - also play FatSat for "short" prepscans that do not have a PASL call --> then play out FatSat for ALL slices, irrespective of FatSat Mode
	if (  
	     !(  lKernelMode == KERNEL_IMAGE )
		  ||  (l_LongEntry_SmsFactor == 1 &&  REOInfo.getPATAccelerationFactorPE() == 1 )
		  ||  ( alPrepScanCounter[lSlice]  > m_lInitialDummyScans + m_lSliceAccPrepScans + m_lPATPrepScans - 1  )
		 )
	  {
cout << " should get FatSat " << endl;
		if (pMrProt->preparationPulses().fatSatMode() == SEQ::FAT_SAT_STRONG  || pMrProt->preparationPulses().fatSatMode() == SEQ::FAT_SAT_WEAK  )
		{
	cout << " getting FatSat " << endl;		
			if(! CSatFat.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice])   )
			{
			  textTR("CSatFat.run failed.");
			  return SEQU_ERROR; 
			}
			if(! SpoilGrad.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice]) )
			{
			  textTR("SpoilGrad.run failed.");
			  return SEQU_ERROR; 
			}
		}
	}
*/

  #endif  //PASL
  
  //---------------------------------------------------------------------------
  // Debug
  //---------------------------------------------------------------------------
  #ifdef SHOW_LOOP_STRUCTURE
    if(mIsDebugBitMask(DEBUG_SEQLOOP|DEBUG_LANDMARK) && lKernelMode!=KERNEL_CHECK)
    { 
      static char tLoopText[512];
	  static char tTemp[256];
      sprintf(tLoopText,"- SBBEPIKernel (%ld echos) with %s",REOInfo.getEchoTrainLength(),ptAddRTEB);
	  sprintf(tTemp,"- slice position: %7f",asSLC[lSlice].getSliceShift());
	  strcat(tLoopText,tTemp);
      mySeqLoop.ShowEvent (tLoopText);
    }
  #endif
  
  
  // benpos SMS-CAIPI -->
	long lOrigSlcIndex = lSlice; // needed to reset after EPI readout
	
	if ( (l_LongEntry_SmsFactor > 1.0) && (mySeqLoop.getInnerSliceNumber() == pMrProt->sliceSeries().size()/l_LongEntry_SmsFactor) )  // this second condition can only be fulfilled if this a multi-banded kernel call 
    {
		//recalculate the slice-index 
		lSlice = getSmsRecalculatedSliceIndex( pMrProt, &asSLC[lSlice],  lSlice,  l_LongEntry_SmsFactor);

		EPIKernel.getReadOutAddress()->Mdh.setCslc(lSlice);
    }

	// <-- benpos SMS-CAIPI

  
  	
  //---------------------------------------------------------------------------
  // execute EPI Kernel
  //---------------------------------------------------------------------------
  if(! EPIKernel.run(pMrProt, pSeqLim, pSeqExpo, &asSLC[lSlice])) /*! EGA-07 !*/
  {
      TRACE_PUT1_NLS(TC_INFO, TF_SEQ, "%s: EPIKernel.run failed.",ptModule,EPIKernel.getNLSStatus());
      return EPIKernel.getNLSStatus();
  }

  // ggg
  gSpiralRO(pMrProt,pSeqLim,pSeqExpo,lKernelMode,lSlice,lPartition,lLine,0,paramLongInterleaves,1,asSLC[lSlice],0, 0 + 0);
  // ggg for SKOPE: Put delay into wipmem
  pMrProt->wipMemBlock().alFree[50]=500;

  //---------------------------------------------------------------------------
  // execute PASL fill times 
  //---------------------------------------------------------------------------

  #ifdef PASL   



  // benpos --> add VASO_BOLD_filltime: only for VASO_BOLD mode, after the last slice of the first "volume per inversion"
	if  ( 
		 	 lKernelMode != KERNEL_PREPARE 
		  && pMrProt->Asl().Mode() == SEQ::ASL_CUSTOM2 
		  && lOrigSlcIndex == pMrProt->sliceSeries().size()/l_LongEntry_SmsFactor - 1
		  && ( EPIKernel.getReadOutAddress()->Mdh.getCrep() ) % 2 == 0 
		)
	{
		fRTEBInit(&asSLC[lSlice].m_sROT_MATRIX);
		fRTEI( l_VasoBoldFilltime   , 0, 0, 0, 0, 0, 0, 0); 
		fRTEBFinish();
	}
	// benpos <--




    if( bEnableFirstPrepScanAsM0Scan )
    {
      // disable readout for further PrepScans after M0 reference scan
      if(lKernelMode==KERNEL_PREPARE)
      {
        if( alPaslPrepScanCounter[lSlice] == 1 )
        {
          fRTSetReadoutEnable(0);
        }
	 }
    }


	if (  
		 (  lKernelMode != KERNEL_PREPARE  && (EPIKernel.getReadOutAddress()->Mdh.getCrep()+ 1 ) % l_LongEntry_VolumesPerInversion == 0   ) // always label for imaging scans
      || (  lKernelMode == KERNEL_PREPARE  && l_LongEntry_SmsFactor == 1 &&  REOInfo.getPATAccelerationFactorPE() == 1 )													// always label for dummies if iPat and SliceAcc are off 
	  || (  lKernelMode == KERNEL_PREPARE  && alPrepScanCounter[lSlice]  > m_lInitialDummyScans + m_lSliceAccPrepScans + m_lPATPrepScans - 1  )
    )
	{
		if( lSlice >= mySeqLoop.getInnerSliceNumber() - 1 ) //benpos ASL
		{

		 NLS_STATUS lStatus = SEQU__NORMAL;
		  //lStatus = fSBBFillTimeRun( lTRfillPASL * pMrProt->sliceSeries().size() );
		  
		 long filltime = lTRfillPASL * mySeqLoop.getInnerSliceNumber(); 
		 
		 if (lKernelMode == KERNEL_PREPARE) 
		 {
			filltime += mySeqLoop.getlKernelScanTime()/l_LongEntry_VolumesPerInversion * (l_LongEntry_VolumesPerInversion - 1) * mySeqLoop.getInnerSliceNumber();
		 }

 	 		 
		 lStatus = fSBBFillTimeRun( filltime ); //benpos ASL
		  if (lStatus != SEQU__NORMAL) 	return (lStatus);
	  }
  }
  #endif
 
    
  
  //---------------------------------------------------------------------------
  // execute mandatory fill time
  //---------------------------------------------------------------------------
  NLS_STATUS lStatus = fSBBCoolTimeRun ( pMrProt, lCoolingPause() );
  if (!lStatus)
  {
      TRACE_PUT1_NLS(TC_INFO, TF_SEQ, "%s: fSBBFillTimeRun failed.",ptModule,lStatus);
      return lStatus;
  }

  // benpos SMS-CAIPI -->
  if (l_LongEntry_SmsFactor > 1.0)
  {
      lSlice = lOrigSlcIndex; //set back to original slice index
  }
  // <-- benpos SMS-CAIPI

  //---------------------------------------------------------------------------
  // increase prep-scan-counter for current slice
  //---------------------------------------------------------------------------
  #ifdef SUPPORT_iPAT_a_ep2d
  {
    if (lKernelMode==KERNEL_PREPARE)
    {
        alPrepScanCounter[lSlice]++;

        if (alPrepScanCounter[lSlice]==mySeqLoop.getlPreparingScans()  )  
        {
            alPrepScanCounter[lSlice]=0;
        }
		//benpos also increment the counter needed for FLEET
		lPrepScanExcitationCounter++;
    }
  }
  #endif

  //---------------------------------------------------------------------------
  // disable read-outs again (may have been switched on for PAT reference scans)
  //---------------------------------------------------------------------------
  #ifdef SUPPORT_iPAT_a_ep2d
  {
    if (lKernelMode==KERNEL_PREPARE)
    {
        fRTSetReadoutEnable(0);
    }
  }
  #endif

  

  //---------------------------------------------------------------------------
  // finished
  //---------------------------------------------------------------------------
  mTRend;
  return SEQU_NORMAL;
}

#ifdef SUPPORT_PACE
NLS_STATUS fSEQReceive(SeqLim* pSeqLim, SeqExpo* pSeqExpo, SEQData& rSEQData)
{
    if( !mySeqLoop.receive(pSeqLim,pSeqExpo,rSEQData) )
    {
        TRACE_PUT2(TC_ALWAYS,TF_SEQ,"Error at %s(%d)",__FILE__,__LINE__);
        return mySeqLoop.getNLSStatus();
    }
    return SEQU_NORMAL;
}
#endif

//-------------------------------------------------------------------------------------
// Helper function for cooling pause calculation
//
// - EPIKernel has to be prepared
// - REOInfo has to be prepared
//
// For some Magnetoms, an elevated gradient performance is used during the
// EPI readout that requires the introduction of an additional pause afterwards
// in order to stay within hardware limitations. This is realized by introducing 
// a mandatory fill time after each kernel. This time is calculated based on the 
// duration of the echo train. 
//
// Exception #1: for the EP2D_SE variant no additional fill time is neccessary
// (no gradient activity between excitation and refocussing)
//
// Exception #2: for the EP2D_DIFF variant, an additional dynamic component is
// required which is calculated based on the gradient load of the diffusion
// module. Note: In this case,  the same pause duration has to be used in 
// fPrepDIFFPlugInForEPIKernel!
//-------------------------------------------------------------------------------------
#ifdef EP2D_DIFF
long lCoolingPause (bool bDynamic)
#else
long lCoolingPause (bool /* bDynamic */)
#endif
{
#ifdef EP2D_SE
    return 0;
#endif
   // Do not use protocol parameters here! Otherwise, one might run
   // into inconsistencies if fSeqPrep is called within a TryProt
   // context!
   const double dIncTRRatio  = theGPA.getTRIncrementEPI();    // TR increment [% / 100]
   long         lEchoSpacing = EPIKernel.getEchoSpacing();
   long         lEchoNumber  = REOInfo.getEchoTrainLength();

   // Static contribution: increment TR by a given ratio of the EPI readout train
   long         lTRIncrement = fSDSRoundDownGRT ( dIncTRRatio * lEchoNumber * lEchoSpacing);

#ifdef EP2D_DIFF
   // Dynamic contribution: increment TR by the absolute value provided by the diffusion SBB
   if ( bDynamic )
   {
    lTRIncrement += fSDSRoundDownGRT ( Diff->getTRIncrement () );
   }
#endif

   // TRACE_PUT1 (TC_INFO, TF_SEQ, "TR increases by %ldus", lTRIncrement);

   return lTRIncrement;
}



// benpos caipi --> helper function for slice index calculation
long getSmsRecalculatedSliceIndex(MrProt* pMrProt, sSLICE_POS* pSLC, long lOrigSlcIndex, long lSmsFactor)
{

	long lNewSliceIndex = 0;

		if    (pMrProt->sliceSeries().isAscending() )
		{
			lNewSliceIndex = lOrigSlcIndex; 
		}
		else if (pMrProt->sliceSeries().isDescending()  ) 
		{
			lNewSliceIndex = lOrigSlcIndex + (lSmsFactor-1) * pMrProt->sliceSeries().size() / lSmsFactor ;
		}
		else if (pMrProt->sliceSeries().isInterleaved()  ) 
		{	
			//interleaved always really means *interleaved ascending* order
			//for an odd  number of slices it goes 0-2-4-6... 1-3-5     --> for an odd number slices the first "interleave" will have one more slice than the seccond
			//for an even number of slices it goes 1-3-5-7... 0-2-4-6   --> for an even number of slices, the first and second "interleave" of slices will have an equal number 
			// --> we can take the ceil value to work out the number of slices for the first "interleave" (true for both the of the full set of slices and slice-accelerated set of slices)
			long lNSlcFirstInterleaveFull = ceil( pMrProt->sliceSeries().size()/2  );
			long lNSlcFirstInterleaveSMS = ceil (pMrProt->sliceSeries().size() / lSmsFactor / 2 );
			long lAnatomicalSliceNo = pSLC->getSliceIndex();

			// rejig the slice index if  the "anatomical range" and "slice index range" are exceeded
			if ( ( lAnatomicalSliceNo >= pMrProt->sliceSeries().size()/lSmsFactor )  && (lOrigSlcIndex >= lNSlcFirstInterleaveSMS ) ) 
			{
				lNewSliceIndex = lOrigSlcIndex - lNSlcFirstInterleaveSMS + lNSlcFirstInterleaveFull  ;
			}
			else
			{
				lNewSliceIndex = lOrigSlcIndex;
			}
		
		}
		else 
		{ 
			cout << " unknown sliceSerriesMode --> I guess we are screwed " << endl; 
		} 

 
		//cout << endl;
		//cout << " Original Slice Index     " << lOrigSlcIndex  << endl;
		//cout << " Modified ChronSliceIndex " << lNewSliceIndex << endl;
		//cout << " Modified AnatSliceIndex  " << asSLC[lNewSliceIndex].getSliceIndex() << endl;

	
		return lNewSliceIndex;
 }    
// <-- benpos caipi