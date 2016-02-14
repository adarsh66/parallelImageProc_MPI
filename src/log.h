/* ------------------------------------------------------------
 * References
 * The code below was adapted from the original on the 
 * website LearnCTheHardWay, author: Zed Shaw
 * Original code can be found at
 * http://c.learncodethehardway.org/book/ex20.html
 * ------------------------------------------------------------
 * debug (msg) 		: if you want to print out from every rank
 * log_it(rank, msg) 	: if you want to print out from rank 0 only.
 */

#ifndef __dbg_h__
#define __dbg_h__

#include <stdio.h>
#include <errno.h>
#include <string.h>

#ifdef NDEBUG
#define debug(M, ...)
#else
#define debug(M, ...) fprintf(stderr, "[DEBUG] (%s:%d): " M, __FILE__, __LINE__, ##__VA_ARGS__)
#endif

#define log_info(M, ...) fprintf(stderr, "[INFO] " M, ##__VA_ARGS__)

#define log_it(rank, M, ...) if(rank==0) {log_info(M,##__VA_ARGS__);}

#endif
