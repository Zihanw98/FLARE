/**
 * @file user_config.h
 * @brief Problem configuration and callback interfaces.
 */
#pragma once

#ifndef HD_NDIM
#define HD_NDIM 1
#endif

// Keep settings already supplied by a problem's user_config.h.
#ifndef HD_NCOMP
#define HD_NCOMP 0
#endif

#ifndef HD_NUSER
#define HD_NUSER 0
#endif
