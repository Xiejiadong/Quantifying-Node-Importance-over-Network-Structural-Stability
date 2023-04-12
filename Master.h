#pragma once

#include <omp.h>
#include <unordered_map>
#include <random>
#include "Partition.h"
#include "Maintenance.h"

class Master
{
public:

	Partition* partition;

	Maintenance* maintenance;

	Master(Partition* partition, Maintenance* maintenance);

	void Anchoring(int b);

	void dynamicMaintain();
};