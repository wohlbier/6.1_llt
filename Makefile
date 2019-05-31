HDRS = algebra.hh types.hh
SRCS = main.cc
EMU_OBJS = $(subst .cc,.emu.o,$(SRCS))

#EMU_PATH = /local/devel/packages/emu-18.11-cplus
EMU_PATH = /local/devel/packages/emu-19.02
#EMU_PATH = /home/jgwohlbier/devel/packages/emu-19.02
EMU_CXX = $(EMU_PATH)/bin/emu-cc
EMU_SIM = $(EMU_PATH)/bin/emusim.x

EMU_SIM_ARGS =
#EMU_SIM_ARGS += --short_trace
#EMU_SIM_ARGS += --memory_trace
EMU_SIM_ARGS += --log2_nodelet_memory 34 --core_clk_mhz 175 --log2_nodelets_per_node 2 --gcs_per_nodelet 3 --ddr_speed 2

EMU_PROFILE = $(EMU_PATH)/bin/emusim_profile

CPPFLAGS =
CPPFLAGS += -D__PROFILE__
LDFLAGS = -lemu_c_utils

EXE  = llt
EMU_EXE = $(EXE).mwx
#INPUT = tri-8-10-0.bin
#INPUT = tri-8-11-1.bin
#INPUT = tri-8-12-3.bin
#INPUT = tri-8-13-5.bin
#INPUT = tri-8-14-7.bin
#INPUT = tri-16-24-0.bin
#INPUT = tri-16-25-2.bin
#INPUT = tri-16-26-4.bin
#INPUT = tri-16-28-8.bin
#INPUT = tri-19-34-9.bin
#INPUT = tri-32-78-63.bin
#INPUT = tri-64-191-184.bin
#INPUT = tri-128-388-379.bin
#INPUT = tri-256-934-994.bin
#INPUT = tri-512-1737-1582.bin
#INPUT = tri-1021-3606-3190.bin
#INPUT = tri-1024-3631-3223.bin
#INPUT = tri-2048-7802-8116.bin
INPUT = triangle_count_data_ca-HepTh-9877-25973-28339.bin

$(EMU_EXE) : $(EMU_OBJS)
	$(EMU_CXX) -o $(EMU_EXE) $(EMU_OBJS) $(LDFLAGS)

run : $(EMU_EXE)
	$(EMU_SIM) $(EMU_SIM_ARGS) $(EMU_EXE) ./tris/$(INPUT)

profile : $(EMU_EXE)
	$(EMU_PROFILE) profile $(EMU_SIM_ARGS) -- $(EMU_EXE) ./tris/$(INPUT)

convert : convert.cc
	$(CXX) -o convert convert.cc

%.emu.o: %.cc $(HDRS)
	$(EMU_CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY : clean

clean :
	-$(RM) *~ $(OBJS) $(EMU_OBJS) $(EXE) $(EMU_EXE) *.cdc *.hdd *.vsf
	-$(RM) -r profile $(EXE).txt
	-$(RM) -r convert
