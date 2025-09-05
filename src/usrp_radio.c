#include "../include/usrp_radio.h"
#include <stdio.h>
#include <string.h>
#include <uhd/types/tune_request.h>
#include <uhd/usrp/usrp.h>

/* This function connects to a USRP b210. And captures a block of samples at
 * sample rate (fs) and centre frequency (fc) */
int usrpIQcapture(float complex *rxSig, const float fs, const float fc,
                  float gain, const unsigned int length) {

  // Create USRP handle
  uhd_usrp_handle usrp;
  if (uhd_usrp_make(&usrp, "") != UHD_ERROR_NONE) {
    fprintf(stderr, "Failed to create USRP\n");
    return EXIT_FAILURE;
  }

  // Set Master clock rate
  uhd_usrp_set_master_clock_rate(usrp, 61.44e6, 0);

  // Set sample rate
  if (uhd_usrp_set_rx_rate(usrp, fs, 0) != UHD_ERROR_NONE) {
    fprintf(stderr, "Failed to set sample rate\n");
    return EXIT_FAILURE;
  }

  // Set centre frequency
  uhd_tune_request_t tune_request = {
      .target_freq = fc,
      .rf_freq_policy = UHD_TUNE_REQUEST_POLICY_AUTO,
      .dsp_freq_policy = UHD_TUNE_REQUEST_POLICY_AUTO};
  uhd_tune_result_t tune_res;
  uhd_usrp_set_rx_freq(usrp, &tune_request, 0, &tune_res);
  if (tune_res.actual_rf_freq == 0.0) {
    fprintf(stderr, "Failed to set centre frequency!\n");
    return EXIT_FAILURE;
  }

  // Set gain
  if (uhd_usrp_set_rx_gain(usrp, gain, 0, "") != UHD_ERROR_NONE) {
    fprintf(stderr, "Failed to set gain\n");
    return EXIT_FAILURE;
  }

  // Create sample streamer
  uhd_stream_args_t stream_args;
  memset(&stream_args, 0, sizeof(stream_args));
  stream_args.cpu_format = "fc32";
  stream_args.otw_format = "sc16";
  stream_args.channel_list = (size_t[]){0};
  stream_args.n_channels = 1;
  stream_args.args = "";

  uhd_rx_streamer_handle rx_stream;
  uhd_rx_streamer_make(&rx_stream);
  uhd_usrp_get_rx_stream(usrp, &stream_args, rx_stream);

  // Stream command
  uhd_stream_cmd_t stream_cmd = {.stream_mode =
                                     UHD_STREAM_MODE_NUM_SAMPS_AND_DONE,
                                 .num_samps = length,
                                 .stream_now = true};

  uhd_rx_streamer_issue_stream_cmd(rx_stream, &stream_cmd);

  // Point to address of rx signal.
  void *buffs[] = {rxSig};

  // Receive samples.
  size_t samples_received;
  uhd_rx_metadata_handle md;
  uhd_rx_metadata_make(&md);
  uhd_rx_streamer_recv(rx_stream, buffs, length, &md, 5.0, false,
                       &samples_received);
  if (samples_received != (size_t)length) {
    printf("Unable to capture the desired no. of samples!\n");
    printf("Captured %d out of %d samples\n", (unsigned int)samples_received,
           length);
    return EXIT_FAILURE;
  }

  // Free resources.
  uhd_usrp_free(&usrp);
  uhd_rx_streamer_free(&rx_stream);
  uhd_rx_metadata_free(&md);

  return EXIT_SUCCESS;
}
