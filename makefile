RM := rm -rf

ITS_SOURCES += ./its/ITS.cpp ./its/main_ITS.cpp 

RMCRAG_SOURCES += ./rmcrag/rmcrag.cpp 

ITS_OBJECTS += $(ITS_SOURCES:.cpp=.o)

RMCRAG_OBJECTS += $(RMCRAG_SOURCES:.cpp=.o)

ITS_EXECUTABLE = ./its/its
RMCRAG_EXECUTABLE = ./rmcrag/rmcrag

.PHONY: all ITS RMCRAG

all: ITS RMCRAG

ITS: $(ITS_EXECUTABLE)

RMCRAG: $(RMCRAG_EXECUTABLE)

$(ITS_EXECUTABLE): $(ITS_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: C++ Linker'
	g++  -o "ITS" $(ITS_OBJS)
	@echo 'Finished building target: $@'
	@echo ' '

$(RMCRAG_EXECUTABLE): $(RMCRAG_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: C++ Linker'
	g++  -o "rmcrag" $(RMCRAG_OBJS)
	@echo 'Finished building target: $@'
	@echo ' '

.cpp.o:
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I./its -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

clean:
	-$(RM) $(ITS_OBJS)$(RMCRAG_OBJS)$(EXECUTABLES)
	-@echo ' '

