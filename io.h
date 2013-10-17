/*
 * NEEMP - io.h
 *
 * by Tomas Racek (tom@krab1k.net)
 * 2013
 *
 * */

#ifndef __IO_H__
#define __IO_H__

#include "subset.h"

void load_molecules(void);
void load_charges(void);
void load_parameters(struct kappa_data * const ss);

void output_charges(const struct subset * const ss);
void output_charges_stats(const struct subset * const ss);
void output_parameters(const struct subset * const ss);

#endif /* __IO_H__ */
