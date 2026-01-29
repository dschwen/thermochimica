#include "thermochimica/context/PhaseConstraints.hpp"
#include <cmath>
#include <algorithm>

namespace Thermochimica {

void PhaseConstraints::allocate(int nSolnPhasesSys, int nConPhasesSys) {
    solnPhaseConstraints.resize(nSolnPhasesSys);
    condPhaseConstraints.resize(nConPhasesSys);

    // Initialize all constraints to unconstrained
    for (auto& c : solnPhaseConstraints) {
        c.mode = PhaseConstraintMode::None;
        c.targetFraction = 0.0;
        c.lagrangeMultiplier = 0.0;
        c.currentFraction = 0.0;
    }
    for (auto& c : condPhaseConstraints) {
        c.mode = PhaseConstraintMode::None;
        c.targetFraction = 0.0;
        c.lagrangeMultiplier = 0.0;
        c.currentFraction = 0.0;
    }
}

void PhaseConstraints::reset() {
    // Reset solver state but keep constraint definitions
    penaltyParameter = 1.0;
    currentOuterIteration = 0;

    // Reset Lagrange multipliers and current fractions
    for (auto& c : solnPhaseConstraints) {
        c.lagrangeMultiplier = 0.0;
        c.currentFraction = 0.0;
    }
    for (auto& c : condPhaseConstraints) {
        c.lagrangeMultiplier = 0.0;
        c.currentFraction = 0.0;
    }
}

void PhaseConstraints::clear() {
    // Clear all constraints
    for (auto& c : solnPhaseConstraints) {
        c.mode = PhaseConstraintMode::None;
        c.targetFraction = 0.0;
        c.lagrangeMultiplier = 0.0;
        c.currentFraction = 0.0;
    }
    for (auto& c : condPhaseConstraints) {
        c.mode = PhaseConstraintMode::None;
        c.targetFraction = 0.0;
        c.lagrangeMultiplier = 0.0;
        c.currentFraction = 0.0;
    }

    // Reset parameters
    penaltyParameter = 1.0;
    penaltyGrowthRate = 10.0;
    constraintTolerance = 1e-4;
    maxOuterIterations = 20;
    currentOuterIteration = 0;
}

bool PhaseConstraints::hasActiveConstraints() const {
    for (const auto& c : solnPhaseConstraints) {
        if (c.mode != PhaseConstraintMode::None) {
            return true;
        }
    }
    for (const auto& c : condPhaseConstraints) {
        if (c.mode != PhaseConstraintMode::None) {
            return true;
        }
    }
    return false;
}

int PhaseConstraints::getNumActiveConstraints() const {
    int count = 0;
    for (const auto& c : solnPhaseConstraints) {
        if (c.mode == PhaseConstraintMode::Fixed) {
            ++count;
        }
    }
    for (const auto& c : condPhaseConstraints) {
        if (c.mode == PhaseConstraintMode::Fixed) {
            ++count;
        }
    }
    return count;
}

bool PhaseConstraints::areConstraintsSatisfied() const {
    for (const auto& c : solnPhaseConstraints) {
        if (c.mode == PhaseConstraintMode::Fixed) {
            double violation = std::abs(c.currentFraction - c.targetFraction);
            if (violation > constraintTolerance) {
                return false;
            }
        }
    }
    for (const auto& c : condPhaseConstraints) {
        if (c.mode == PhaseConstraintMode::Fixed) {
            double violation = std::abs(c.currentFraction - c.targetFraction);
            if (violation > constraintTolerance) {
                return false;
            }
        }
    }
    return true;
}

double PhaseConstraints::getTotalConstraintViolation() const {
    double totalViolation = 0.0;
    for (const auto& c : solnPhaseConstraints) {
        if (c.mode == PhaseConstraintMode::Fixed) {
            double violation = c.currentFraction - c.targetFraction;
            totalViolation += violation * violation;
        }
    }
    for (const auto& c : condPhaseConstraints) {
        if (c.mode == PhaseConstraintMode::Fixed) {
            double violation = c.currentFraction - c.targetFraction;
            totalViolation += violation * violation;
        }
    }
    return totalViolation;
}

} // namespace Thermochimica
