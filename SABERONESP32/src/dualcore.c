/*
 * dualcore.c
 *
 *  Created on: 2019年3月7日
 */

#include "dualcore.h"
#include "stdio.h"
#include "stdint.h"
#include "cpucycles.h"
//#include "freertos/portmacro.h"
//#include "freertos/task.h"
//#include "freertos/projdefs.h"

#define EXT_TASK_PRIORITY  10

void dualcoreInitSemaphore() {
	if (semaphore_main == NULL) {
		semaphore_main = xSemaphoreCreateCounting(SEMAPHORE_MAX_COUNT, 0);
	}
	if (semaphore_ext == NULL) {
		semaphore_ext = xSemaphoreCreateCounting(SEMAPHORE_MAX_COUNT, 0);
	}

}

void dualcoreCreateTask(TaskFunction_t task, void * const params) {
	xTaskCreatePinnedToCore(task, "ext_task", 20 * 1024, params, EXT_TASK_PRIORITY, NULL, CORE_EXT);
}

void dualcoreTakeSemaphore(SemaphoreHandle_t s) {
#if DUALCORE_DEBUG_SEMAPHORE
	//dualcoreDebug(s == semaphore_main ? "main" : " ext", "Mul A[i]*sk ready");
	printf("Take ");
	printf(s == semaphore_main ? "main" : "ext");
	printf("\n");
#endif
	xSemaphoreTake(s, portMAX_DELAY);

}

void dualcoreGiveSemaphore(SemaphoreHandle_t s) {
#if DUALCORE_DEBUG_SEMAPHORE
	//dualcoreDebug(s == semaphore_main ? "main" : " ext", "Mul A[i]*sk ready");
	printf("Give ");
	printf(s == semaphore_main ? "main" : "ext");
	printf("\n");
#endif
	xSemaphoreGive(s);

}

void dualcoreDeleteTask() {
	vTaskDelete(NULL);
}

void dualcoreDebug(char* taskName, char* msg) {
	//uxTaskGetStackHighWaterMark(xTask)
	printf("Core %u [%4s] [%12u] : %s\n", xPortGetCoreID(), taskName, cpucycles(), msg);
}

