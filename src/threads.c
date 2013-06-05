/** Types for Multi-threading */

/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Copyright (C) 2012 Genome Research Ltd.                                  *
 *                                                                           *        
 *  Author: Hannes Ponstingl (hp3@sanger.ac.uk)                              *
 *                                                                           *
 *  This file is part of SMALT.                                              *
 *                                                                           *
 *  SMALT is free software: you can redistribute it and/or modify it under   *
 *  the terms of the GNU General Public License as published by the Free     *
 *  Software Foundation, either version 3 of the License, or (at your        *
 *  option) any later version.                                               *
 *                                                                           *
 *  This program is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU         *
 *  General Public License for more details.                                 *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License along  *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.          *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************/
#include <stdlib.h>
#include <semaphore.h>
#include <pthread.h>

#include "threads.h"

enum TRHEAD_CONST {
  THREAD_NUM_MAX = 64,  /**< Maximum number of threads to be spawned */
  BUFFARG_NUM_FAC_DEFAULT  = 8, /**< Number of buffered thread arguments/results
				 * as a multiplicative factor. I.e. for n>0 threads there are
				 * n*BUFFARG_NUM_FAC buffered arguments */
  THREAD_BUFF_NUM = 3,  /**< Number of buffers for arguments that are passed between
			 * threads. */
  THREAD_TASK_NUM = 4,   /**< Number of tasks available */
};

enum THREADS_STATUS_FLAGS {
  THRFLG_INIT = 0x01,  /**< Threads were initialised */
  THRFLG_PROC = 0x02,
  THRFLG_INPUT = 0x04,
  THRFLG_OUTPUT = 0x08,
  THRFLG_ARGBUF = 0x10, /**< Argument buffers have been set up */
  THRFLG_SETUP = 0x20,   /**< Threads have been set up */
  THRFLG_STARTED = 0x40
};

enum THREAD_ARGUMENT_FLAGS {
  THRARG_INDEPT = 0x01,
  THRARG_SIGNED = 0x02,
};

enum THREAD_BUFFER_TYPES {
  THRBUFTYP_EMPTY,
  THRBUFTYP_LOADED,
  THRBUFTYP_PROCESSED,
  THRBUFTYP_UNKNOWN
};

typedef uint8_t BOOL_t;
typedef uint8_t THRSTATUSFLG_t;
typedef uint8_t THRARGFLG_t;
typedef void *(PTHREAD_WRAPPER_FUNC)(void *);

typedef struct _BUFFARG { /**< Arguments passed between threads */
#ifdef threads_debug
  short threadno;
  uint64_t readno;
#endif
  short argno;
  void *thisp;
  struct _BUFFARG *nextp;
} BUFFARG;

typedef struct _ARGBUFF { /**< Buffer for arguments to be passed between threads */
  int nThreadsPushing;  /**< Number of threads pushing on this buffer,
			 * 0 signals termination */
  BUFFARG *firstp;
  BUFFARG *lastp;
  sem_t sema;            /* semaphore, keeps track of how many arguments are in buffer */
  pthread_mutex_t mutex; /* lock access to shared data */
#ifdef threads_debug
  uint8_t buftyp;
#endif
} ARGBUFF;

typedef struct _THREADTASK {
  THRSTATUSFLG_t status;
  short n_threads;      /**< Number of threads for this task (can be 0) */
  THREAD_INITF *initf;   /**< Initialisation function */
  const void *initargp;  /**< Additional argument to pass to intialisation function */
  THREAD_PROCF *procf;   /**< Processing function */
  THREAD_CLEANF *cleanf; /**< Function for cleaning up */
  THREAD_CHECKF *checkf; /**< Checks buffer before fetching */
  THREAD_CMPF *cmpf;     /**< Compare buffer arguments */
  size_t argsz;         /**< Size of the thread argument */
  uint8_t fromx;        /**< Index of buffer from which arguments are pulled */
  uint8_t tox;          /**< Index of buffer to which argments are pushed */
} THREADTASK;

typedef struct _THREADARG {
  short threadno;   /**< Thread number (< 0 means not run as a separate thread */
  THRARGFLG_t flags; /**< Combination of THREAD_ARGUMENT_FLAGS */
  uint8_t task;     /**< one of THREAD_TASKS */
  void *p;          /**< points to thread-specific data (SmaltMapArgs or SmaltInputArgs) */
  BUFFARG *buflstp;  /**< Start of linked list for internal buffering (e.g. sorted ouptput) */
  ErrMsg *errmsgp;  /**< Thread specific error messages */
  int exit_code; /**< error code with which wrapper exits */
} THREADARG;

static struct _Threads {
  uint8_t status;       /**< One of THREADS_STATUS_FLAGS */
  short n_threads;      /**< Number of threads (0 if run single-threaded) */
  THREADTASK tasks[THREAD_TASK_NUM];
  pthread_t *threadp;   /**< Thread structures, array of size n_threads */
  short n_targ;         /**< Size of array targp */
  THREADARG *targp;     /**< Thread arguments array of size n_targ */
  short n_buffargs;     /**< Number of bufferend arguments */
  BUFFARG *buffargp;    /**< Buffered arguments, array of size n_buffargs */
  void *memp;           /**< Memory allocated for arguments */
  ARGBUFF buff[THREAD_BUFF_NUM];/**< Buffers for arguments that are passed between threads.
  				 * buffp[0]: empty arguments, bufp[1]: loaded arguments.
  				 * buffp[2]: processed arguments. */
} Threads;

#ifdef threads_debug
static pthread_mutex_t mutex_stdout = PTHREAD_MUTEX_INITIALIZER;
#endif

#ifdef threads_debug
static int getBufArgNum(const ARGBUFF *fifop)
{
  int n = 0;
  BUFFARG *hp;
  for (hp = fifop->firstp; NULL != hp; hp = hp->nextp, n++);
  return n;
}

static const char *getBufTyp(const ARGBUFF *fifop)
{
  switch (fifop->buftyp) {
  case THRBUFTYP_EMPTY:
    return "empty";
  case THRBUFTYP_LOADED:
    return "loaded";
  case THRBUFTYP_PROCESSED:
    return "processed";
  case THRBUFTYP_UNKNOWN:
  default:
    return "unknown";
  }
}

static const char *getTaskTyp(const THREADARG *p)
{
  switch(p->task) {
  case THRTASK_ARGBUF:
    return "argument";
  case THRTASK_INPUT:
    return "input";
  case THRTASK_PROC:
    return "processing";
  case THRTASK_OUTPUT:
    return "output";
  default:
    return "unknown";
  }
}
#endif

/*****************************************************************************
 *********************** Methods of private type ARGBUFF *********************
 *****************************************************************************/

static int pushARGBUFF(ARGBUFF *fifop, BUFFARG *argp)
/**< Put argument in buffer. argp == NULL initialises. */
{
  int errcode = ERRCODE_SUCCESS;
#ifdef threads_debug
  int ns = 0;
  int na = 0;
#endif

  if (argp == NULL) { /* initialisation signal */
    pthread_mutex_lock(&fifop->mutex);
    fifop->firstp = fifop->lastp = NULL;
    fifop->nThreadsPushing = 0;
    pthread_mutex_unlock(&fifop->mutex);
#ifdef threads_debug
  pthread_mutex_lock(&mutex_stdout);
  fprintf(stderr, "#THREAD_DEBUG:initialised '%s' buffer ...\n", 
	  getBufTyp(fifop));
  pthread_mutex_unlock(&mutex_stdout);
#endif

  } else {
    pthread_mutex_lock(&fifop->mutex);
    if (fifop->firstp == NULL) {
      if (fifop->lastp != NULL)
	errcode = ERRCODE_ASSERT;
      fifop->firstp = fifop->lastp = argp;
      argp->nextp = NULL;
    } else {
      if (fifop->lastp == NULL)
	errcode = ERRCODE_ASSERT; 
      fifop->lastp = fifop->lastp->nextp = argp;
      argp->nextp = NULL;       
    }
#ifdef threads_debug
    na = getBufArgNum(fifop);
#endif
    pthread_mutex_unlock(&fifop->mutex);
    sem_post(&fifop->sema);
#ifdef threads_debug
  sem_getvalue(&fifop->sema, &ns);
  pthread_mutex_lock(&mutex_stdout);
  printf("THREAD_DEBUG:thread[%hi]: pushed to '%s' buffer: sema = %i, nThreadsPushing = %i , n_args = %i, errcode = %i ...\n", 
	 argp->threadno, getBufTyp(fifop), ns, fifop->nThreadsPushing, na, errcode);
  pthread_mutex_unlock(&mutex_stdout);
#endif

  }

  return errcode;
}

static int pullARGBUFF(BUFFARG **argp, 
		       ARGBUFF *fifop)
/**< Fetch argument from buffer. Returns termination code 
 *  ERRCODE_PTHRTERMSIG if no argument in buffer and fifop->nThreadsPushing == 0.
 */
{
  int errcode = ERRCODE_SUCCESS;
#ifdef threads_debug
  int ns = 0;
  int na = 0;
#endif
#ifdef threads_debug
  sem_getvalue(&fifop->sema, &ns);
  pthread_mutex_lock(&mutex_stdout);
  printf("THREAD_DEBUG:pull from '%s' buffer: waiting (sema = %i, nThreadsPushing = %i) ...\n", 
	 getBufTyp(fifop), ns, fifop->nThreadsPushing);
  pthread_mutex_unlock(&mutex_stdout);
#endif

  sem_wait(&fifop->sema);

#ifdef threads_debug
  sem_getvalue(&fifop->sema, &ns);
  pthread_mutex_lock(&mutex_stdout);
  printf("THREAD_DEBUG:pull from '%s' buffer: finished waiting (sema = %i, nThreadsPushing = %i) ...\n", 
	 getBufTyp(fifop), ns, fifop->nThreadsPushing);
  pthread_mutex_unlock(&mutex_stdout);
#endif

  pthread_mutex_lock(&fifop->mutex);
  *argp = fifop->firstp;
  if ((*argp) == NULL) {
    if (fifop->nThreadsPushing <= 0)
      errcode = ERRCODE_PTHRTERMSIG;
    else if (fifop->lastp != NULL)
      errcode = ERRCODE_ASSERT;
  } else {
    if (fifop->firstp == fifop->lastp) {
      if ((*argp)->nextp != NULL)
	errcode = ERRCODE_ASSERT;
      fifop->firstp = fifop->lastp = NULL;
    } else {
      fifop->firstp = (*argp)->nextp;
      (*argp)->nextp = NULL;
    }
  }
#ifdef threads_debug
  na = getBufArgNum(fifop);
#endif
  pthread_mutex_unlock(&fifop->mutex);
  if (errcode == ERRCODE_PTHRTERMSIG)
    sem_post(&fifop->sema); /* release other threads waiting to pull */

#ifdef threads_debug
  sem_getvalue(&fifop->sema, &ns);
  pthread_mutex_lock(&mutex_stdout);
  printf("THREAD_DEBUG:finished pulling from '%s' buffer: sema = %i, nThreadsPushing = %i, n_args = %i, errcode = %i ...\n", 
	 getBufTyp(fifop), ns, fifop->nThreadsPushing, na, errcode);
  pthread_mutex_unlock(&mutex_stdout);
#endif

  return errcode;
}

static int initARGBUFF(ARGBUFF *fifop)
/**< Initialise buffer */
{
  int errcode;
  sem_init(&fifop->sema, 0, 0);
  pthread_mutex_init(&fifop->mutex, NULL);
  errcode = pushARGBUFF(fifop, NULL);
  return errcode;
}

static int signOnARGBUFF(ARGBUFF *fifop)
/**< Increment the counter for threads pushing onto the buffer and
 * return the number */
{
  int n;
  pthread_mutex_lock(&fifop->mutex);
  n = ++(fifop->nThreadsPushing);
  pthread_mutex_unlock(&fifop->mutex);
#ifdef threads_debug
  pthread_mutex_lock(&mutex_stdout);
  printf("THREAD_DEBUG:signed on buffer '%s' (%i) ...\n", getBufTyp(fifop), n);
  pthread_mutex_unlock(&mutex_stdout);
#endif
  return n;
}

static int signOffARGBUFF(ARGBUFF *fifop)
/**< Decrement the counter for threads pushing onto the buffer.
 * Return the number of threads still pushing onto the buffer.
 * If there is no thread pushing, wake threads waiting to pull from the buffer.
 */
{
  int n;
  pthread_mutex_lock(&fifop->mutex);
  n = --(fifop->nThreadsPushing);
  pthread_mutex_unlock(&fifop->mutex);
  if (n <= 0) {
     sem_post(&fifop->sema);
  }
  
#ifdef threads_debug
  pthread_mutex_lock(&mutex_stdout);
  printf("THREAD_DEBUG:signed off buffer '%s' (%i) ...\n", getBufTyp(fifop), n);
  pthread_mutex_unlock(&mutex_stdout);
#endif
  
  return n;
}

static void termARGBUFF(ARGBUFF *fifop)
/**< Set termination signal for buffer */
{
#ifdef threads_debug
  pthread_mutex_lock(&mutex_stdout);
  printf("threads_debug:termARGBUFF() called ...\n");
  pthread_mutex_unlock(&mutex_stdout);
#endif
  pthread_mutex_lock(&fifop->mutex);
  fifop->nThreadsPushing = 0;
  pthread_mutex_unlock(&fifop->mutex);
  sem_post(&fifop->sema);
}

/****************************************************************************
 ******************** Methods of Private Type THREADARG *********************
 ****************************************************************************/

static void pushTHREADARGInternalBuffer(THREADARG *thargp, BUFFARG *argp, 
				       THREAD_CMPF *cmpf)
{
  BUFFARG *hp = thargp->buflstp;
#ifdef threads_debug
  BUFFARG *p;
  int i = 0;
  if (NULL != cmpf) {
    pthread_mutex_lock(&mutex_stdout);
  
    printf("THREADS_DEBUG: +++++ pushTHREADARGInternalBuffer: threadno = %hi, readno = %llu +++++\n",
	   thargp->threadno, (unsigned long long) argp->readno); 
    for (p = thargp->buflstp, i = 0; (p); p = p->nextp, i++)
      printf("THREADS_DEBUG: pushTHREADARGInternalBuffer: [%i] %llu\n",
	     i, (unsigned long long) p->readno);
    pthread_mutex_unlock(&mutex_stdout);
  }
#endif
  if (NULL == cmpf || NULL == hp || 
      cmpf(argp->thisp, hp->thisp) <= 0) {
    argp->nextp = hp;
    thargp->buflstp = argp;
  } else {
    BUFFARG *lp = hp->nextp;
    for(; (lp) && cmpf(argp->thisp, lp->thisp) > 0; lp = lp->nextp)
      hp = lp;
    argp->nextp = lp;
    hp->nextp = argp;
  }
#ifdef threads_debug
  if (NULL != cmpf) {
    pthread_mutex_lock(&mutex_stdout);
    printf("THREADS_DEBUG: + pushTHREADARGInternalBuffer: on exit: threadno = %hi +\n", thargp->threadno);
    for (p = thargp->buflstp, i = 0; (p); p = p->nextp, i++)
      printf("THREADS_DEBUG: pushTHREADARGInternalBuffer: readno[%i] = %llu\n",
	     i, (unsigned long long) p->readno);
    printf("THREADS_DEBUG: ++ pushTHREADARGInternalBuffer: exit ++\n");
    pthread_mutex_unlock(&mutex_stdout);
  }
#endif
}

static BUFFARG *pullTHREADARGInternalBuffer(THREADARG *thargp, THREAD_CHECKF *checkf, void *tdatap)
{
  BUFFARG *argp = thargp->buflstp;
#ifdef threads_debug
  int chkrv = 0;
  BUFFARG *p;
  int i = 0;
  if (NULL != checkf) {
    pthread_mutex_lock(&mutex_stdout);
    printf("THREADS_DEBUG: ----- pullTHREADARGInternalBuffer threadno = %hi-----\n", thargp->threadno); 
    for (p = thargp->buflstp, i = 0; (p); p = p->nextp, i++)
      printf("THREADS_DEBUG: pullTHREADARGInternalBuffer: [%i] %llu\n",
	     i, (unsigned long long) p->readno);
    pthread_mutex_unlock(&mutex_stdout);
  }
#endif

  if (argp != NULL) {

    if (NULL == checkf ||
#ifdef threads_debug
	(chkrv = checkf(tdatap, argp->thisp))) {
#else
      (checkf(tdatap, argp->thisp))) {
#endif
      thargp->buflstp = argp->nextp;
      argp->nextp = NULL;
    } else {
      argp = NULL;
    }
  }

#ifdef threads_debug
  if (NULL != checkf) {
    pthread_mutex_lock(&mutex_stdout);
    printf("THREADS_DEBUG: pullTHREADARGInternalBuffer: checkf = %i\n",
	   chkrv);
    printf("THREADS_DEBUG: - pullTHREADARGInternalBuffer: on exit: threadno = %hi -\n", thargp->threadno);
    for (p = thargp->buflstp, i = 0; (p); p = p->nextp, i++)
      printf("THREADS_DEBUG: pullTHREADARGInternalBuffer: readno[%i] = %llu\n",
	     i, (unsigned long long) p->readno);
    printf("THREADS_DEBUG: -- pullTHREADARGInternalBuffer: exit --\n");
    pthread_mutex_unlock(&mutex_stdout);
  }
#endif

  return argp;
}
/****************************************************************************
 *************** PTHREAD_WRAPPER_FUNC function run by threads ***************
 ****************************************************************************/

static void *tprocf(void *p) 
{
  int errcode = ERRCODE_SUCCESS;
  BUFFARG *argp;
  THREADARG *thargp = (THREADARG *) p;
  THREADTASK *taskp = Threads.tasks + thargp->task;
  ARGBUFF *argbf_fromp = Threads.buff + taskp->fromx;
  ARGBUFF *argbf_top = Threads.buff + taskp->tox;

#ifdef threads_debug
  pthread_mutex_lock(&mutex_stdout);
  fprintf(stderr, "# THREADS_DEBUG:thead[%hi]: tprocf(): task '%s' started\n", 
	 thargp->threadno, getTaskTyp(thargp));
  pthread_mutex_unlock(&mutex_stdout);
#endif

  if ( !(thargp->flags & THRARG_SIGNED)) {
#ifdef threads_debug
    pthread_mutex_lock(&mutex_stdout);
    fprintf(stderr, 
	    "# THREADS_DEBUG:thead[%hi]:tprocf(): task '%s': sign on buffer typ '%s'\n", 
	   thargp->threadno, getTaskTyp(thargp), getBufTyp(Threads.buff + taskp->tox));
    pthread_mutex_unlock(&mutex_stdout);
#endif
    signOnARGBUFF(argbf_top);
    thargp->flags |= THRARG_SIGNED;
  }

  do {

    if ((errcode = pullARGBUFF(&argp, argbf_fromp)))
      break;

#ifdef threads_debug
    pthread_mutex_lock(&mutex_stdout);
    fprintf(stderr, "# THREAD_DEBUG:thead[%hi]:tprocf(): pulled argument %i\n", 
	   thargp->threadno, argp->argno);
    pthread_mutex_unlock(&mutex_stdout);
#endif
    if (argp == NULL) {
      errcode = ERRCODE_ASSERT;
      break;
    } 
    
#ifdef threads_debug
    pthread_mutex_lock(&mutex_stdout);
    fprintf(stderr, 
	    "# THREAD_DEBUG:thead[%hi]:tprocf('%s', %hu): push to internal buffer" \
	   "argument %i, read = %llu, cmpf %c= NULL\n", 
	   thargp->threadno, getTaskTyp(thargp), 
	   (unsigned short) thargp->task, 
	   argp->argno, 
	   (unsigned long long) argp->readno, (NULL == taskp->cmpf)? '=':'!');
    pthread_mutex_unlock(&mutex_stdout);
#endif

    pushTHREADARGInternalBuffer(thargp, argp, taskp->cmpf);

    while ((ERRCODE_SUCCESS == errcode) && 
	   (argp = pullTHREADARGInternalBuffer(thargp, taskp->checkf, thargp->p))) {
      errcode = taskp->procf(thargp->errmsgp, 
#ifdef threads_debug
			     &argp->readno,
#endif
			     thargp->p, argp->thisp);

    /* push to loaded arguments */
#ifdef threads_debug
      pthread_mutex_lock(&mutex_stdout);
      fprintf(stderr,
	      "# THREADS_DEBUG:thead[%hi]:tprocf push read onto buffer %hi\n", 
	   thargp->threadno, taskp->tox);
      pthread_mutex_unlock(&mutex_stdout);
#endif
      if (!(errcode))
	errcode = pushARGBUFF(argbf_top, argp);
    }
 
  } while ((ERRCODE_SUCCESS == errcode)
	   && thargp->threadno >= 0);

  if ((errcode)) { 
#ifdef threads_debug
    pthread_mutex_lock(&mutex_stdout);
    fprintf(stderr,
	    "# THREADS_DEBUG:thead[%hi]:sign off from buffer %hi.\n", 
	   thargp->threadno, taskp->tox);
    pthread_mutex_unlock(&mutex_stdout);
#endif

    signOffARGBUFF(argbf_top);
    thargp->flags &= ~THRARG_SIGNED;
  }
#ifdef threads_debug
  pthread_mutex_lock(&mutex_stdout);
  fprintf(stderr,
	  "# THREADS_DEBUG:thead[%hi]:tprocf(): task '%s' finished.\n", 
	 thargp->threadno, getTaskTyp(thargp));
  pthread_mutex_unlock(&mutex_stdout);
#endif
  thargp->exit_code = errcode;

  return p;
}

/****************************************************************************
 ***************************** Public Methods *******************************
 ****************************************************************************/

void threadsInit(void)
{
  short i;
  Threads.status = 0;
  Threads.threadp = NULL;
  Threads.targp = NULL;
  Threads.buffargp = NULL;
  for (i=0; i<THREAD_BUFF_NUM; i++) {
#ifdef threads_debug
    if (i == 0)
      Threads.buff[i].buftyp = THRBUFTYP_EMPTY;
    else if (i == 1)
      Threads.buff[i].buftyp = THRBUFTYP_LOADED;
    else if (i == 2)
      Threads.buff[i].buftyp = THRBUFTYP_PROCESSED;
    else
      Threads.buff[i].buftyp = THRBUFTYP_UNKNOWN;     
#endif
    initARGBUFF(Threads.buff + i);
  }
  Threads.status |= THRFLG_INIT;
}

int threadsSetTask(uint8_t task_typ, 
		   short n_threads, 
		   THREAD_INITF *initf, 
		   const void *initargp,
		   THREAD_PROCF *procf, 
		   THREAD_CLEANF *cleanf,
		   THREAD_CHECKF *checkf,
		   THREAD_CMPF *cmpf,
		   size_t argsz)
{
  int errcode = ERRCODE_SUCCESS;
  THREADTASK *taskp = NULL;

  if (task_typ >= THREAD_TASK_NUM)
    return ERRCODE_FAILURE;

  if (!(Threads.status & THRFLG_INIT))
    return ERRCODE_ASSERT;

  taskp = Threads.tasks + task_typ;
  taskp->n_threads = n_threads;
  taskp->initf = initf;
  taskp->initargp = initargp;
  taskp->procf = procf;
  taskp->cleanf = cleanf;
  taskp->checkf = checkf;
  taskp->cmpf = cmpf;
  taskp->argsz = argsz;
  taskp->status = THRFLG_INIT;
  taskp->fromx = taskp->tox = 0;

  switch (task_typ) {
  case THRTASK_INPUT:
    taskp->fromx = THRBUFTYP_EMPTY;
    taskp->tox = THRBUFTYP_LOADED;
    Threads.status |= THRFLG_INPUT;
    break;
  case THRTASK_ARGBUF:
    taskp->n_threads = 0;
    taskp->procf = NULL;
    Threads.status |= THRFLG_ARGBUF;
    break;
  case THRTASK_PROC:
    taskp->fromx = THRBUFTYP_LOADED;
    taskp->tox = THRBUFTYP_PROCESSED;
    Threads.status |= THRFLG_PROC;
    break;
  case THRTASK_OUTPUT:
    taskp->fromx = THRBUFTYP_PROCESSED;
    taskp->tox = THRBUFTYP_EMPTY;
    Threads.status |= THRFLG_OUTPUT;
#ifdef threads_debug
    fprintf(stderr,
	    "#THREADS_DEBUG: threadsSetTask [%hu]: cmpf %c= NULL, checkf %c= NULL\n",
	   (unsigned short) task_typ, 
	   (NULL ==  taskp->cmpf)? '=':'!', 
	   (NULL ==  taskp->checkf)? '=':'!');
#endif
    break;
  default:
      errcode = ERRCODE_FAILURE;
    break;
  }

  return errcode;
}

int threadsSetUp(int n_buffarg_factor)
{
  int errcode = ERRCODE_SUCCESS;
  uint8_t testflg = THRFLG_INIT | THRFLG_PROC | THRFLG_INPUT | THRFLG_OUTPUT | THRFLG_ARGBUF;
  short i, nta, nth, ntr, n_threads, n_targ;
  size_t memsz;
  void *hp;

  if (n_buffarg_factor <= 0)
    n_buffarg_factor = BUFFARG_NUM_FAC_DEFAULT;
  if ((Threads.status & testflg) != testflg || 
      ((testflg & THRFLG_SETUP)))
    return ERRCODE_ASSERT;

  Threads.n_threads = 0;
  Threads.n_targ = 0;
  memsz = 0;
  for (nta=0; nta<THREAD_TASK_NUM; nta++) {
    if (nta != THRTASK_ARGBUF) {
      THREADTASK *taskp = Threads.tasks + nta;
      if (taskp->n_threads < 1) {
	Threads.n_targ += 1;
	memsz += taskp->argsz;	
      } else {
	Threads.n_threads += taskp->n_threads;
	Threads.n_targ += taskp->n_threads;
	memsz += taskp->argsz * taskp->n_threads;
      }
    }
  }

  Threads.n_buffargs = (Threads.n_threads > 0)? Threads.n_threads * n_buffarg_factor: 1;
  memsz += Threads.n_buffargs * Threads.tasks[THRTASK_ARGBUF].argsz;

  if (Threads.n_threads > 0) {
    ECALLOCP(Threads.n_threads, Threads.threadp);
    if (NULL == Threads.threadp)
      errcode = ERRCODE_NOMEM;
  } else {
    Threads.threadp = NULL;
  }
    
  ECALLOCP(Threads.n_targ, Threads.targp);  
  ECALLOCP(Threads.n_buffargs, Threads.buffargp);
  Threads.memp = EMALLOC0(memsz);
    
  if (NULL == Threads.targp || 
      NULL == Threads.buffargp ||
      NULL == Threads.memp)
    errcode = ERRCODE_NOMEM;

  n_threads = 0;
  n_targ = 0;
  hp = Threads.memp;
  for (nta=0, nth=0, ntr=0; nta<THREAD_TASK_NUM && !(errcode); nta++) {
    if (nta != THRTASK_ARGBUF) {
      THREADTASK *taskp = Threads.tasks + nta;
      if (taskp->n_threads < 1) {
	n_targ++;
      } else {
	n_threads += taskp->n_threads;
	n_targ += taskp->n_threads;
      }
      for (; ntr<n_targ && !(errcode); ntr++) {
	THREADARG *targp = Threads.targp + ntr;
	targp->task = nta;
	targp->p = hp;
	targp->buflstp = NULL;
	hp += taskp->argsz;
	if (taskp->n_threads < 1) {
	  targp->threadno = -1;
	  targp->flags = 0;
	} else {
	  targp->threadno = nth++;
	  targp->flags = THRARG_INDEPT;
	}
	if (NULL == ERRMSG_CREATE(targp->errmsgp))
	  errcode = ERRCODE_NOMEM;
	else
	  errcode = (*Threads.tasks[nta].initf)(targp->p, taskp->initargp, targp->threadno);
      }
    }
  }
  if (!(errcode) && 
      (Threads.n_threads != n_threads ||
       Threads.n_targ != n_targ))
    errcode = ERRCODE_ASSERT;

  for (i=0;i<Threads.n_buffargs && (!errcode); i++) {
    BUFFARG *p = Threads.buffargp + i;
    p->thisp = hp;
    hp += Threads.tasks[THRTASK_ARGBUF].argsz;
    p->argno = i;
    p->nextp = NULL;
#ifdef threads_debug
    p->threadno = 0;
    p->readno = 0;
#endif
    errcode = (*Threads.tasks[THRTASK_ARGBUF].initf)(p->thisp, 
						      Threads.tasks[THRTASK_ARGBUF].initargp, 
						      i);
    if (ERRCODE_SUCCESS == errcode) {
      errcode = pushARGBUFF(&Threads.buff[THRBUFTYP_EMPTY], p);
    }
  }

  if (ERRCODE_SUCCESS == errcode)
    Threads.status |= THRFLG_SETUP;

  return errcode;
}

void threadsCleanup(void)
{
  if (Threads.status & THRFLG_SETUP) {
    short i;
    for (i=0;i<Threads.n_buffargs; i++) {
      BUFFARG *bargp = Threads.buffargp + i;
      (*Threads.tasks[THRTASK_ARGBUF].cleanf)(NULL, bargp->thisp);
    }
    for (i=0; i<Threads.n_targ; i++) {
      THREADARG *targp = Threads.targp + i;
      (*Threads.tasks[targp->task].cleanf)(targp->errmsgp, targp->p);
      ERRMSG_END(targp->errmsgp)
    }
    free(Threads.memp);

    free(Threads.buffargp);
    free(Threads.targp);
    free(Threads.threadp);
  }
  Threads.status = 0;
}

int threadsStart(void)
{
  int rv = ERRCODE_SUCCESS;
  short ntr;

  if (!(Threads.status & THRFLG_SETUP))
    return ERRCODE_ASSERT;

  /* start up thread for input */ 
  for (ntr=0; ntr < Threads.n_targ; ntr++) {
    THREADARG *targp = Threads.targp + ntr;
    if (targp->threadno < 0)
      continue;
    if ((rv = pthread_create(Threads.threadp + targp->threadno, NULL, 
			     tprocf, targp)))
      ERRMSGNO(targp->errmsgp, ERRCODE_PTHREAD);
  }
  if (!rv) 
    Threads.status |= THRFLG_STARTED;

  return rv;
}

void threadsStop(void)
{
  int nth;
    /* wait for threads to finish */
  for (nth=0; nth < Threads.n_threads; nth++) {
#ifdef threads_debug
    pthread_mutex_lock(&mutex_stdout);
    fprintf(stderr,
	    "# THREAD_DEBUG: main_thread: waiting for thread %hi to finish ...\n", nth);
    pthread_mutex_unlock(&mutex_stdout);
#endif
    pthread_join(Threads.threadp[nth], NULL);
  }

  return;
}

int threadsRun(void)
{
  int errcode = ERRCODE_SUCCESS;
  BOOL_t has_nonthread = 1;

  threadsStart();

  
  while (ERRCODE_SUCCESS == errcode && (has_nonthread)) {
    short ntr;
    has_nonthread = 0;
    for(ntr=0; ntr<Threads.n_targ && (ERRCODE_SUCCESS == errcode); ntr++) {
      THREADARG *targp = Threads.targp + ntr;
      if (targp->threadno < 0) {
	has_nonthread = 1;
	tprocf(targp);
	errcode = targp->exit_code;
      }
    }
  }
   
  threadsStop();

  return (ERRCODE_PTHRTERMSIG == errcode)? ERRCODE_SUCCESS: errcode;
}

void *threadsGetMem(uint8_t task_typ)
{
  void *p = NULL;

  if (Threads.status & THRFLG_SETUP) {
    if (task_typ == THRTASK_ARGBUF) {
      p = Threads.buffargp[0].thisp;
    } else {
      short ntr;
      for (ntr=0; ntr < Threads.n_targ && p == NULL; ntr++) {
	if (Threads.targp[ntr].task == task_typ) 
	  p = Threads.targp[ntr].p;
      }
    }
  }
  
  return p;
}