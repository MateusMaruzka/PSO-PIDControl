/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  * @attention
  *
  * <h2><center>&copy; Copyright (c) 2020 STMicroelectronics.
  * All rights reserved.</center></h2>
  *
  * This software component is licensed by ST under BSD 3-Clause license,
  * the "License"; You may not use this file except in compliance with the
  * License. You may obtain a copy of the License at:
  *                        opensource.org/licenses/BSD-3-Clause
  *
  ******************************************************************************
  */
/* USER CODE END Header */

/* Includes ------------------------------------------------------------------*/
#include "main.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */

#include "arm_math.h"

#include <stdlib.h>


/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */

#define NDIM 2
#define SWARM 10

#define ITER_MAX 50
/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/

/* USER CODE BEGIN PV */

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
/* USER CODE BEGIN PFP */

float genRandomNumber(void); // ???

void fillVector(float *vector, unsigned int len);

float esfera(float x, float y);

unsigned int argmin(float *vector, unsigned int len);

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

float genRandomNumber(void){
	//static int i = 0;
	//return i++;
	return (float)rand() / (float)RAND_MAX - 0.5;
}


void fillVector(float *vector, unsigned int len){
	
	for(int i = 0; i < len; i++){
		
		*(vector+i) = genRandomNumber();
		
	}
	
}



float esfera(float x, float y){
	
	return x*x + y*y;
	
}


unsigned int argmin(float *vector, unsigned int len){

	unsigned int idx = 0;
	float min = 0;
	
	for(unsigned int i = 0; i < len; i++){
		if (*(vector+i) < min){
			idx = i;
		}
	}
	
	return idx;
}
	


void updateVelocity(arm_matrix_instance_f32 X, 
										arm_matrix_instance_f32 V, 
										arm_matrix_instance_f32 pBest,
										float32_t *gBest,
										uint32_t ndim){
											
	// v = v*w + c1*r1(pbest-x) + c2*r2(gbest - x)
											
}
/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */

int main(void)
{
  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  /* USER CODE BEGIN 2 */	

	// Swarm`s position and velicity
	float32_t x[SWARM][NDIM], 
						v[SWARM][NDIM], 
						fitness[SWARM],
						fitBest[SWARM],
						pBest[SWARM][NDIM],
						aux[SWARM][NDIM];
	
	// Initialize the matrices
	fillVector((float32_t *)x, NDIM*SWARM);
	fillVector((float32_t *)v, NDIM*SWARM);
	
	float32_t c1 = 4, c2 = 4;
	
	// pBest recebe a posicao do enxame
	arm_copy_f32((const float32_t *)x, (float32_t *)pBest, NDIM*SWARM);

	arm_matrix_instance_f32 X, V, FITNESS, PBEST, AUX, FITBEST;

	arm_mat_init_f32(&X, SWARM, NDIM, (float32_t *)x);
	arm_mat_init_f32(&V, SWARM, NDIM, (float32_t *)v);
	
	arm_mat_init_f32(&FITNESS, SWARM, 1, fitness);
	arm_mat_init_f32(&FITBEST, SWARM, 1, fitBest);
	
	arm_mat_init_f32(&PBEST, SWARM, NDIM, (float32_t *)pBest);
	arm_mat_init_f32(&AUX, SWARM,NDIM, (float32_t *)aux);
	
	
	for(unsigned int i = 0; i < SWARM; i++){
		arm_power_f32((float32_t *)(x+i), NDIM, fitness+i);
	}
	
	arm_mat_scale_f32(&FITNESS, 1.0, &FITBEST);


	// Função Esfera implementada de outra maneira

	float32_t fitgBest;
	uint32_t idx_gb = 0;
	arm_min_f32(fitBest, SWARM, &fitgBest, &idx_gb);
	
	for(int i = 0; i < ITER_MAX; i++){
		
		
		// Calculando a velocidade de cada partícula
		arm_mat_scale_f32(&V, 0.5, &V);
		
		arm_mat_sub_f32(&PBEST, &X, &AUX);
		
		arm_mat_scale_f32(&AUX, 0.5, &AUX);
		
		arm_mat_add_f32(&V, &AUX, &V);
	
		for(unsigned int j = 0; j < SWARM; j++){
			arm_sub_f32(pBest+idx_gb*NDIM, x+j, aux+j, NDIM);
			//arm_power_f32((float32_t *)(x+i), NDIM, aux+i);
		}
		
		arm_mat_scale_f32(&AUX, 0.5, &AUX);

		arm_mat_add_f32(&V, &AUX, &V);
		
		// Fim calculo velocidade
		
		// Incrementa posicao
		arm_mat_add_f32(&X, &V, &X);	
		
		
		// Função Esfera implementada de outra maneira
		for(int k = 0; k < SWARM; k++){
			arm_power_f32((float32_t *)(x+k), NDIM, fitness+k);
		}
		
		// Atualiza pbest
		for(unsigned int l = 0; l < SWARM; l++){
			
			if(fitness[l] <  fitBest[l] || (int)fitBest[l] == -1){
				
				fitBest[l] = fitness[l];
				arm_copy_f32(pBest + l*NDIM, x+l*NDIM, NDIM);
				
			}
			
		}
		// Encontra gbest
		arm_min_f32(fitBest, SWARM, &fitgBest, &idx_gb);
	
	
	}
	
	
  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */

  while (1) {
    /* USER CODE END WHILE */		
    /* USER CODE BEGIN 3 */
  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Initializes the CPU, AHB and APB busses clocks 
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI;
  RCC_OscInitStruct.HSIState = RCC_HSI_ON;
  RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_NONE;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }
  /** Initializes the CPU, AHB and APB busses clocks 
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_HSI;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_0) != HAL_OK)
  {
    Error_Handler();
  }
}

/* USER CODE BEGIN 4 */

/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */

  /* USER CODE END Error_Handler_Debug */
}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{ 
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     tex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */

/************************ (C) COPYRIGHT STMicroelectronics *****END OF FILE****/
