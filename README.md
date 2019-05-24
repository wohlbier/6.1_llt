Count triangles using L L^T.

Significantly change the parallelism from 6_llt. Spawn thread to each nodelet
and have each nodelet responsible for rows it owns.
