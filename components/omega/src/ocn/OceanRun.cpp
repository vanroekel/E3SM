//===-- ocn/OceanRun.cpp - Run Ocean Model ----------------------*- C++ -*-===//
//
// The ocnRun method advances the model forward from CurrTime until the
// EndAlarm rings.
//
//===----------------------------------------------------------------------===//

#include "Config.h"
#include "Forcing.h"
#include "IOStream.h"
#include "OceanDriver.h"
#include "OceanState.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include <filesystem>

namespace OMEGA {

static bool forcingInputFileAvailable() {
   Error Err;

   Config *OmegaConfig = Config::getOmegaConfig();
   Config StreamsConfig("IOStreams");
   Err += OmegaConfig->get(StreamsConfig);
   if (Err.isFail())
      return false;

   Config ForcingConfig("Forcing");
   Err += StreamsConfig.get(ForcingConfig);
   if (Err.isFail())
      return false;

   bool UsePointerFile = false;
   Err += ForcingConfig.get("UsePointerFile", UsePointerFile);
   if (Err.isFail() || UsePointerFile)
      return true;

   std::string Filename;
   Err += ForcingConfig.get("Filename", Filename);
   if (Err.isFail() || Filename.find("$") != std::string::npos)
      return true;

   return std::filesystem::exists(Filename);
}

int ocnRun(TimeInstant &CurrTime ///< [inout] current sim time
) {

   // error code
   I4 Err = 0;

   // fetch default OceanState and TimeStepper
   OceanState *DefOceanState   = OceanState::getDefault();
   TimeStepper *DefTimeStepper = TimeStepper::getDefault();
   Forcing *DefForcing         = Forcing::getDefault();

   // EndAlarm must be set before calling ocnRun
   OMEGA_REQUIRE(DefTimeStepper->hasEndAlarm(), "ocnRun: no EndAlarm");

   // get simulation time and other time info
   Clock *OmegaClock     = DefTimeStepper->getClock();
   Alarm *EndAlarm       = DefTimeStepper->getEndAlarm();
   TimeInterval TimeStep = DefTimeStepper->getTimeStep();
   TimeInstant SimTime   = OmegaClock->getCurrentTime();

   // Get Simulation metadata field for later updates
   std::shared_ptr<Field> SimInfo = Field::get(SimMeta);

   // time loop, integrate until EndAlarm or error encountered
   I8 IStep = 0;
   while (Err == 0 && !(EndAlarm->isRinging())) {

      // track step count
      ++IStep;

      // Refresh optional file-based forcing fields if the Forcing stream
      // exists and is scheduled to read at this model time.
      if (forcingInputFileAvailable()) {
         Metadata ForcingReqMeta;
         Error ForcingReadErr =
             IOStream::read("Forcing", OmegaClock, ForcingReqMeta);
         if (ForcingReadErr.isFail()) {
            if (ForcingReadErr.Msg.find("Stream Forcing not found") ==
                std::string::npos) {
               CHECK_ERROR(ForcingReadErr,
                           "Errors encountered reading Forcing during run");
               ABORT_ERROR(
                   "Error updating forcing variables from input stream");
            }
         }
      }

      // call forcing routines, anything needed pre-timestep
      DefForcing->computeAll();

      // do forward time step
      // first call to doStep can sometimes take very long
      // we want to time it separately and disable child timers
      // for that timer
      if (IStep == 1) {
         Pacer::start("Stepper:firstDoStep", 1);
         Pacer::disableTiming();
         DefTimeStepper->doStep(DefOceanState, SimTime);
         Pacer::enableTiming();
         Pacer::stop("Stepper:firstDoStep", 1);
      } else {
         Pacer::start("Stepper:doStep", 1);
         DefTimeStepper->doStep(DefOceanState, SimTime);
         Pacer::stop("Stepper:doStep", 1);
      }

      // write restart file/output, anything needed post-timestep

      IOStream::writeAll(OmegaClock);

      LOG_INFO("ocnRun: Time step {} complete, clock time: {}", IStep,
               SimTime.getString(4, 4, "-"));
   }

   return Err;

} // end ocnRun

} // end namespace OMEGA
