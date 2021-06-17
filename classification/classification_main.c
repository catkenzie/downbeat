// File name: TestOverlapAdd.c
//
// 27 Sep 2016 .. coding started .. KM
// 24 Jan 2018 .. updated for EECS 452 STM32F7xx .. KM

#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "arm_math.h"
#include "arm_const_structs.h"

#include "hamming.h"
#include "b_fixed.h"

#define FFT_SIZE 1024
#define NUM_SPECBANDS 25
#define NUM_KIT_PIECES 3

void Startup(void);

// 	ADC FIFO call definitions
uint16_t ADC_Fifo_Init(void);
uint16_t ADC_Fifo_Get(int32_t *);

// DAC FIFO call definitions
uint16_t DAC_Fifo_Init(void);
uint16_t DAC_Fifo_Put(int32_t);

void LED_Green_On();

static uint16_t ctr;

static int32_t ADC_input;
static float *ptr1, *ptr2, *fftptr;

// Input buffers, split into halves to allow for overlap
static float first_half[FFT_SIZE];
static float second_half[FFT_SIZE];
static float FFT_buff[2*FFT_SIZE];

// Compressed spectrum data of input
static float specband[NUM_SPECBANDS];

// G factor of input resulting from NNMF
static float gain_factors[NUM_KIT_PIECES];

// Variables for Onset Detection Thresholding
static bool detected_hit[NUM_KIT_PIECES] = {false, false, false};
static int onset_levels[NUM_KIT_PIECES] = {0};
static int onset_frames[NUM_KIT_PIECES] = {0};
static float attack_thresh[NUM_KIT_PIECES] = {2, 3.3, 7.7};
static float decay_thresh[NUM_KIT_PIECES] = {0.002, 0.5, 4.3};

// Array storing number of FFT frequency components in each frequency band
uint16_t HZBANDS[] = {3, 2, 2, 3, 2, 3, 3, 4, 4, 4, 5, 5, 7, 7, 9, 11, 12, 17,
					  21, 25, 30, 42, 58, 81, 152};

// For calculating time (1 frame corresponds to 1 FFT, remember overlap)
int frame = 0;

char onset_time[16];

// REQUIRES: FFT_buff points to the first element in a signal in the Hz-domain
//			 Sum of elements in HZBANDS = FFT_SIZE/2
// MODIFIES: specband
// EFFECT: Computes the band-wise sums of the signal frame's spectrum
void hz_band_summation(float* FFT_buff, float* specband)
{
	float real, imag, fft_power;
	int samp = 0;

	for (int hz_band = 0; hz_band < NUM_SPECBANDS; hz_band++)
	{
		specband[hz_band] = 0; //since specband is reused for each frame
		for (int iter = 0; iter < HZBANDS[hz_band]; iter++)
		{
			real = FFT_buff[samp];
			imag = FFT_buff[samp + 1];
			fft_power = sqrt(real*real + imag*imag);
			specband[hz_band] = specband[hz_band] + fft_power;
			samp += 2;
		}
	}
}

// REQUIRES: specband and b_fixed precalculated
// MODIFIES: gain_factors
// EFFECT: Calculates G matrix (array) for input signal using fixed Basis matrix
void nnmf(float* specband, float* gain_factors)
{
	for (int i = 0; i < NUM_KIT_PIECES; i++)
	{
		gain_factors[i] = 1; //since gain_factors is reused for each frame
	}

	for (int iter = 0; iter < 100; iter++)
	{
		// Used to store intermediate calculations
		float ans_intermed_a[NUM_SPECBANDS] = {0};
		float ans_intermed_b[NUM_KIT_PIECES] = {0};
		float ans_intermed_c[NUM_KIT_PIECES] = {0};

		for (int i = 0; i < NUM_SPECBANDS; i++)
		{
			for (int j = 0; j < NUM_KIT_PIECES; j++)
			{
				ans_intermed_a[i] += (b_fixed[i][j] * gain_factors[j]);
			}
			ans_intermed_a[i] = specband[i] / ans_intermed_a[i];
		}

		for (int i = 0; i < NUM_KIT_PIECES; i++)
		{
			for (int j = 0; j < NUM_SPECBANDS; j++)
			{
				ans_intermed_b[i] += (b_fixed[j][i] * ans_intermed_a[j]);
				ans_intermed_c[i] += b_fixed[j][i];
			}
			gain_factors[i] *= (ans_intermed_b[i] / ans_intermed_c[i]);
		}
	}
}

// REQUIRES: fp a pointer to a file that is open and supports write operations
// MODIFIES: fp
// EFFECT: Determines if an onset is true or false for each piece of the kit
//		   If it is true, it writes the time at which it was hit to a file
void onset_thresholding(float* gain_factors)
{
	for (int piece = 0; piece < NUM_KIT_PIECES; piece++)
	{
		// Was the amplitude strong enough to be an onset?
		if (gain_factors[piece] > attack_thresh[piece])
		{
			detected_hit[piece] = true;
			// Has the onset reached its peak?
			if (gain_factors[piece] > onset_levels[piece])
			{
				onset_frames[piece] = frame;
				onset_levels[piece] = gain_factors[piece];
			}
		}
		if (detected_hit[piece] && frame - onset_frames[piece] == 4)
		{
			// How sharply has the onset decayed over 4 frames?
			if (gain_factors[piece] > decay_thresh[piece])
			{
				// Level still high enough to be considered a true onset
				sprintf(onset_time, "%f", (((onset_frames[piece] - 1) * 512.0) / 44100.0));
//				switch (piece)
//				{
//					case 1: fputs("hi-hat,", fp);
//							break;
//					case 2: fputs("kick,", fp);
//							break;
//					case 3: fputs("snare,", fp);
//							break;
//				}
//				fputs(onset_time, fp);
//				fputs("\n", fp);
			}
			detected_hit[piece] = false;
			onset_levels[piece] = 0;
		}
	}
}

int main()
{
    Startup();   // start what this system needs to have started
    // IS THE ABOVE LINE NECESSARY??

	// initialize FIFOs
	DAC_Fifo_Init(); // currently not using DAC, may be needed for metronome??
	ADC_Fifo_Init();

	LED_Green_On();		// indicate successful FIFO initialization

	// Open/create file for transcription, add first line (header)
//	FILE* fp = fopen("/onsets.csv", "w");
//	fputs("instrument,time(s)\n", fp);

	// Initialize first half of FFT buffer array
	ptr1 = &first_half[0];
	for (ctr = 0; ctr < FFT_SIZE / 2; ctr++)
	{
		while (ADC_Fifo_Get(&ADC_input) == 0);
		*ptr1++ = ADC_input>>16;
		*ptr1++ = 0;
	}

	while (1)
	{
		// Read in second half of FFT buffer array
		ptr2 = &second_half[0];
		for (ctr = 0; ctr < FFT_SIZE / 2; ctr++)
		{
			while (ADC_Fifo_Get(&ADC_input) == 0);
			while(DAC_Fifo_Put((int32_t)gain_factors[0]<<16)==0);
			*ptr2++ = ADC_input>>16;
 			*ptr2++ = 0;
		}

		// Store full input (1024 samples) in FFT_buff
		memcpy(FFT_buff, first_half, FFT_SIZE * sizeof(float));
		memcpy(FFT_buff + FFT_SIZE, second_half, FFT_SIZE * sizeof(float));
		frame++;

		// Window the input signal
		fftptr = &FFT_buff[0];
		for (ctr = 0; ctr < FFT_SIZE; ctr++)
		{
			*fftptr++ *= HAM_1024[ctr];
			fftptr++;
		}

		// Take FFT
		arm_cfft_f32(&arm_cfft_sR_f32_len1024, FFT_buff, 0, 1);

		// Perform frequency band summation
		hz_band_summation(FFT_buff, specband);

		// Get the gain values via NNMF
		nnmf(specband, gain_factors);

		// Onset detection thresholding
		onset_thresholding(gain_factors);

		// Shift second half of FFT buffer array to first half
		memcpy(first_half, second_half, FFT_SIZE * sizeof(float));
	}
}
