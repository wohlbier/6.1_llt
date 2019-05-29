HDRS = algebra.hh types.hh utils.hh
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

EMU_PROFILE = $(EMU_PATH)/bin/emusim_profile

CPPFLAGS =
CPPFLAGS += -D__PROFILE__
LDFLAGS = -lemu_c_utils

EXE  = llt
EMU_EXE = $(EXE).mwx
#INPUT = tri-8-10-0.tsv
#INPUT = tri-8-11-1.tsv
#INPUT = tri-8-12-3.tsv
#INPUT = tri-8-13-5.tsv
#INPUT = tri-8-14-7.tsv
#INPUT = tri-16-24-0.tsv
#INPUT = tri-16-25-2.tsv
#INPUT = tri-16-26-4.tsv
#INPUT = tri-16-28-8.tsv
#INPUT = tri-19-34-9.tsv
#INPUT = tri-32-78-63.tsv
#INPUT = tri-64-191-184.tsv
#INPUT = tri-128-388-379.tsv
#INPUT = tri-256-934-994.tsv
INPUT = tri-512-1737-1582.tsv
#INPUT = tri-1021-3606-3190.tsv
#INPUT = tri-1024-3631-3223.tsv
#INPUT = tri-2048-7802-8116.tsv
#INPUT = triangle_count_data_ca-HepTh-9877-25973-28339.tsv

$(EMU_EXE) : $(EMU_OBJS)
	$(EMU_CXX) -o $(EMU_EXE) $(EMU_OBJS) $(LDFLAGS)

run : $(EMU_EXE)
	$(EMU_SIM) $(EMU_SIM_ARGS) $(EMU_EXE) ./tris/$(INPUT)

profile : $(EMU_EXE)
	CORE_CLK_MHZ=175 \
	$(EMU_PROFILE) profile $(EMU_SIM_ARGS) -- $(EMU_EXE) ./tris/$(INPUT)

%.emu.o: %.cc $(HDRS)
	$(EMU_CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY : clean

clean :
	-$(RM) *~ $(OBJS) $(EMU_OBJS) $(EXE) $(EMU_EXE) *.cdc *.hdd *.vsf
	-$(RM) -r profile $(EXE).txt
