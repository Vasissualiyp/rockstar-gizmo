#ifndef _IO_GIZMO_H_
#define _IO_GIZMO_H_
#ifdef ENABLE_HDF5
#include <stdint.h>
#include "../particle.h"

void load_particles_gizmo(char *filename, struct particle **p, int64_t *num_p);

#endif /* ENABLE_HDF5 */
#endif /* _IO_GIZMO_H_ */
