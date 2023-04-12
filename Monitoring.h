#pragma once

#include <omp.h>
#include <unordered_map>
#include <fstream>
#include "Partition.h"

class Monitoring
{
public:

	Partition* partition;

	Monitoring(Partition* partition);

	void monitor();
};