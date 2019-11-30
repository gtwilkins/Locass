PREFIX = /usr/local/bin/
DBG=-g
INC=-Isrc -Isrc/assembly -Isrc/calibrate -Isrc/commands -Isrc/deadaptor -Isrc/index -Isrc/query -Isrc/shared -Isrc/transform -Isrc/correct
VPATH=src:src/assembly:src/calibrate:src/assembly:src/commands:src/deadaptor:src/index:src/query:src/shared:src/transform:src/correct

SRCS =  \
	locass.cpp \
	assemble.cpp \
	calibrate.cpp \
	correct.cpp \
	index.cpp \
	seed.cpp \
	bitfiles.cpp \
	calibrate_writer.cpp \
	calibrate_structs.cpp \
	correct_amplicon.cpp \
	correct_read.cpp \
	correct_query.cpp \
	correct_structs.cpp \
	deadapter.cpp \
	deamplify.cpp \
	deamplify_structs.cpp \
	error.cpp \
	extend.cpp \
	export_file.cpp \
	filenames.cpp \
	import_file.cpp \
	index_reader.cpp \
	index_structs.cpp \
	index_writer.cpp \
	timer.cpp \
	transform_bwt.cpp \
	transform_structs.cpp \
	transform.cpp \
	transform_binary.cpp \
	parameters.cpp \
	shared_functions \
	leap.cpp \
	local_alignment.cpp \
	locus.cpp \
	locus_calibration.cpp \
	locus_export.cpp \
	locus_extension.cpp \
	locus_fill.cpp \
	locus_gap.cpp \
	locus_leaping.cpp \
	locus_path.cpp \
	locus_pathing.cpp \
	locus_pathing_structs.cpp \
	locus_port.cpp \
	node.cpp \
	node_bridging.cpp \
	node_calibration.cpp \
	node_claim.cpp \
	node_completion.cpp \
	node_coverage.cpp \
	node_cloning.cpp \
	node_export.cpp \
	node_extension.cpp \
	node_filling.cpp \
	node_folding.cpp \
	node_furthest.cpp \
	node_groupings.cpp \
	node_groups.cpp \
	node_islands.cpp \
	node_leaping.cpp \
	node_limits.cpp \
	node_looping.cpp \
	node_navigation.cpp \
	node_pairing.cpp \
	node_pairs.cpp \
	node_path.cpp \
	node_pathing.cpp \
	node_pruning.cpp \
	node_reads.cpp \
	node_reliability.cpp \
	node_seed.cpp \
	node_slicing.cpp \
	node_structs.cpp \
	node_validation.cpp \
	node_verify.cpp \
	path_alleles.cpp \
	path_cross.cpp \
	path_merge.cpp \
	path_reassembly.cpp \
	path_review.cpp \
	path_seed.cpp \
	path_sequence.cpp \
	prune_bubble.cpp \
	query.cpp \
	query_binary.cpp \
	query_extension.cpp \
	query_graph.cpp \
	query_junction.cpp \
	query_state.cpp \
	query_structs.cpp \
	seed_fork.cpp
	

# C++ compiler
CXX = g++
# C++ flags; passed to compiler
CXXFLAGS = -std=c++11
# Linker flags; passed to compiler
LDFLAGS = -std=c++11
# Dependency flags; passed to compiler
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
# Objects directory
OBJDIR = .o
$(shell mkdir -p $(OBJDIR) >/dev/null)
# Dependencies directory
DEPDIR = .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
# Derive objects from sources
OBJS = $(patsubst %,$(OBJDIR)/%.o,$(basename $(SRCS)))
# Derive dependencies from sources
DEPS = $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS)))

# Generic link executable
LINK.o = $(CXX) $(LDFLAGS) $(DBG) -o $@
# Generic compile object
COMPILE.cc = $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(DBG) $(INC) -c -o $@
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

all: locass
	@echo
	@echo 'Compile successful. Type "sudo make install" to complete intall.'

.PHONY: install
install:
	@mv -f locass $(PREFIX)
	@make clean
	@echo
	@echo 'Install successful. Type "locass -h" to see usage.'

.PHONY: clean
clean:
	@$(RM) -r $(OBJDIR) $(DEPDIR)

locass: $(OBJS)
	$(LINK.o) $^

$(OBJDIR)/%.o : %.cpp
$(OBJDIR)/%.o : %.cpp $(DEPDIR)/%.d
	$(COMPILE.cc) $<
	$(POSTCOMPILE)

# Create dependency rule
.PRECIOUS = $(DEPDIR)/%.d
$(DEPDIR)/%.d: ;

# Include dependencies; The '-' ensures no errors
-include $(DEPS)
