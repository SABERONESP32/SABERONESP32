/*
 * dualcore.h
 *
 *  Created on: 2019年3月7日
 *      Author: Mai
 */

#ifndef DUALCORE_H_
#define DUALCORE_H_

#include "freertos/FreeRTOS.h"
#include "freertos/task.h"
#include "freertos/projdefs.h"
#include "freertos/portmacro.h"
#include "freertos/semphr.h"
#include "Arduino.h"



//#define DUALCORE_DEBUG  1
#define DUALCORE_DEBUG_SEMAPHORE  0

#ifdef __cplusplus
extern "C" {
#endif

#define SEMAPHORE_MAX_COUNT 24//8 信号量最大值


//typedef void * QueueHandle_t;
//typedef QueueHandle_t SemaphoreHandle_t;

SemaphoreHandle_t semaphore_main, semaphore_ext ;

#define CORE_ARDUINO  1
#define CORE_EXT      0


typedef struct {
	void * p0;
	void * p1;
	void * p2;
	void * p3;
} PARAMS_4_T;

typedef struct {
	void * p0;
	void * p1;
	void * p2;
	void * p3;
	void * p4;
} PARAMS_5_T;

typedef struct {
	void * p0;
	void * p1;
	void * p2;
	void * p3;
	void * p4;
	void * p5;
} PARAMS_6_T;

typedef struct {
	void * p0;
	void * p1;
	void * p2;
	void * p3;
	void * p4;
	void * p5;
	void * p6;
} PARAMS_7_T;

typedef struct {
	void * p0;
	void * p1;
	void * p2;
	void * p3;
	void * p4;
	void * p5;
	void * p6;
	void * p7;
	void * p8;
} PARAMS_9_T;


void dualcoreInitSemaphore();
//切记，task运行完之后必须执行dualcoreDeleteTask()否则会崩溃
void dualcoreCreateTask(TaskFunction_t task, void * const params);

void dualcoreTakeSemaphore(SemaphoreHandle_t s);
void dualcoreGiveSemaphore(SemaphoreHandle_t s);

//delete task itself
void dualcoreDeleteTask();

//debug message output
void dualcoreDebug(char* taskName, char* msg);

#ifdef __cplusplus
}
#endif



#endif /* DUALCORE_H_ */
