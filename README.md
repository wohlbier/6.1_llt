Count triangles using L L^T.

Significantly change the parallelism from 6_llt. Spawn thread to each nodelet
and have each nodelet responsible for rows it owns.



System notes:
Migration Engine (ME) has 28 queues. 14 incoming and 14 outgoing.
Incoming and outgoing queue to each nodelet Nodelet Queue Manager (NQM) (8 of
each queue type). Remaining 6 queues are connected to the system
interconnect.


Configuration Data and Summary Performance Output (*.cdc)

For the ME queues incoming from the NQMs (ME[ ].FromNQM[ ]):
∗ Number of incoming transactions
∗ Incoming bandwidth utilization
∗ Number of packets being sent to other nodelets on the same node (# to nqms)
∗ Number of packets being sent to the system interconnect
∗ Additional data for packet routing through the ME (non-zero only for Chick system model).
– For the ME queues incoming from the system interconnect (ME[ ].FromSysIC[ ]):
∗ Number of incoming transactions
∗ Incoming bandwidth utilization
∗ Outgoing bandwidth utilization
∗ Additional data for packet routing through the ME (non-zero only for Chick system model).
– For the ME queues outgoing to the NQMs (ME[ ].ToNQM[ ]):
∗ Number of outgoing transactions
∗ Outgoing bandwidth utilization
∗ Additional data for packet routing through the ME (number of requesters to this output queue)
– For the ME queues outgoing to the system interconnect (ME[ ].ToSysIC[ ]):
∗ Number of outgoing transactions
∗ Outgoing bandwidth utilization
∗ Additional data for packet routing through the ME (number of requesters to this output queue)

Verbose Statistics Information (*.vsf)
* Maximum number and time(s) of packets in each queue in the NQM
* Additional instruction based statistics:
– Instructions retired
– Number of threads removed from a GC due to a timeout interrupt
– Number of threads rescheduled (user requested)
– Number of service requests made
– Instructions flushed: instructions that issued but could not be completed due to a re- source being unavailable
– Fences retired: these occur when an instruction is ready to execute, but cannot as the thread is still waiting for acknowledgments to return.
* Additional memory statistics:
– Number of memory operations: number of incoming memory requests to the memory
front end
– Number of memory transactions: number of read/write operations in memory. This number is larger than the number of memory operations as atomics and remote opera- tions count as two transactions (one for the read, one for the write)
– Number of releases: should match the number of threads that quit on the nodelet
* GC-NQM statistics: for these statistics, the values in each GC are summed to get the total
– Count of the number of cycles (stall length[ ]) a thread waited in GC before the NQM could service it
– Number of transactions for each packet length (gc to nqm lengths[ ]) leaving the GCs and going to the NQM
– outgoing xactions, outgoing cycles, outgoing bw utilization: totals and utilization infor- mation for transferring data from the GCs to the NQM
– Similar data for the last two above items is collected for NQM to GC transactions
* NQM Total Stats:
– Number of incoming and outgoing transactions for each packet length
– Summary counters and bandwidth utilization for incoming and outgoing packets

Timed Queue Depths (*.tqd)
Sampled queue depths vs time

Memory Tracing (*.mt)
* TID: thread ID
* SRC: source nodelet
* Type: All memory accesses are grouped into one of three types:
– L: local memory access. Remotes to the local nodelet are given this type. Also, threads that are doing their first memory access post-migration are NOT listed here (i.e., the access that caused the migration is only given an M type).
– R: remote memory operation (remote packet generated)
– M: memory access causing a migration
* DEST: destination nodelet. For Type L accesses, this will match SRC.
* @cycle: the cycle number where the access occurred. This value is not reset when switching from functional to timed mode; the timed mode starts its cycle count by adding 1 from the last cycle in the functional mode.
* TPC(0x): Program counter, in hex, of instruction causing the memory access.
* A(0x): Hex value in the A register when the access occurred
* LocalA(0x): Hex value that is the physical memory location being accessed on the local nodelet. This column may eventually be changed in a future revision to be a global physical address


Approach:
Compute maximum possible bandwidths per nodelet
- NCDRAM
- NQM
- ME

Memory DDR4-1600: Bandwidth = 1.416 GiB/s = 1.521 GB/s
SRIO SystemIC bandwidth     = 2.320 GiB/s = 2.500 GB/s

Total threads: 64 x 3 x 4 x 8 = 6048
Profile shows initial jump to near 5000, but then dials back to ~1000. Why
is there a limit on active threads?
