PREFIX = /usr/local/bin/
INC=-Isrc -Isrc/assembly -Isrc/calibrate -Isrc/commands -Isrc/index -Isrc/query -Isrc/shared -Isrc/transform
VPATH=src:src/assembly:src/calibrate:src/assembly:src/commands:src/index:src/query:src/shared:src/transform

SRCS =  \
	locass.cpp \
	index_writer.cpp \
	index_reader.cpp \
	assemble.cpp \
	index.cpp \
	calibrate.cpp \
	transform_bwt.cpp \
	transform_structs.cpp \
	transform.cpp \
	transform_binary.cpp \
	calibrate_writer.cpp \
	calibrate_structs.cpp \
	timer.cpp \
	parameters.cpp \
	filenames.cpp \
	shared_functions \
	node_looping.cpp \
	locus_export.cpp \
	node_folding.cpp \
	node_validation.cpp \
	locus_pathing_structs.cpp \
	node_export.cpp \
	node_structs.cpp \
	seed.cpp \
	node_completion.cpp \
	node_slicing.cpp \
	node_groupings.cpp \
	node_bridging.cpp \
	node_reliability.cpp \
	locus_fill.cpp \
	node_seed.cpp \
	node_cloning.cpp \
	node_pairing.cpp \
	node_islands.cpp \
	node_pathing.cpp \
	node.cpp \
	extend.cpp \
	node_extension.cpp \
	node_navigation.cpp \
	locus_pathing.cpp \
	export_file.cpp \
	locus_extension.cpp \
	node_reads.cpp \
	locus_calibration.cpp \
	node_furthest.cpp \
	locus.cpp \
	node_coverage.cpp \
	import_file.cpp \
	locus_leaping.cpp \
	node_calibration.cpp \
	path_merge.cpp \
	path_reassembly.cpp \
	path_review.cpp \
	path_sequence.cpp \
	query_state.cpp \
	query_binary.cpp \
	query_extension.cpp \
	query.cpp \
	query_structs.cpp


# C++ compiler
CXX = g++
# C++ flags; passed to compiler
CXXFLAGS = -std=c++11 -g
# Linker flags; passed to compiler
LDFLAGS = -std=c++11 -g
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
LINK.o = $(CXX) $(LDFLAGS) -o $@
# Generic compile object
COMPILE.cc = $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(INC) -c -o $@
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
