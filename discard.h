/*
 * NEEMP - discard.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#ifndef __DISCARD_H__
#define __DISCARD_H__

#include "subset.h"

struct subset *discard_iterative(const struct subset * const full);
struct subset *discard_simple(const struct subset * const initial);

#endif /* __DISCARD_H__ */
