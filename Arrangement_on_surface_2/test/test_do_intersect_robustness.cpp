// Standalone Ray-Ray and Ray-Segment Robustness Test
// Uses minimal kernel - no external dependencies

#include "minimal_kernel.h"
#include <cassert>
#include <vector>
#include <string>

struct TestResult {
    std::string name;
    bool passed;
    std::string message;
};

std::vector<TestResult> results;

void test(const std::string& name, bool condition, const std::string& msg = "") {
    TestResult res;
    res.name = name;
    res.passed = condition;
    res.message = msg;
    results.push_back(res);

    std::cout << (condition ? "✓ PASS" : "✗ FAIL") << ": " << name;
    if (!msg.empty()) std::cout << " (" << msg << ")";
    std::cout << "\n";
}

// ============================================================================
// RAY-RAY TESTS
// ============================================================================

void test_ray_ray_same_source() {
    Ray_2 r1(Point_2(0, 0), Point_2(1, 0));
    Ray_2 r2(Point_2(0, 0), Point_2(0, 1));
    test("Ray-Ray: same source", ray_ray_intersect(r1, r2));
}

void test_ray_ray_non_intersecting() {
    Ray_2 r1(Point_2(0, 0), Point_2(1, 0));
    Ray_2 r2(Point_2(2, 0), Point_2(1, 0));  // Points left
    test("Ray-Ray: opposite directions", !ray_ray_intersect(r1, r2));
}

void test_ray_ray_intersecting() {
    Ray_2 r1(Point_2(0, 0), Point_2(1, 1));
    Ray_2 r2(Point_2(0, 2), Point_2(1, -1));
    test("Ray-Ray: V-shape (should intersect)", ray_ray_intersect(r1, r2));
}

void test_ray_ray_parallel() {
    Ray_2 r1(Point_2(0, 0), Point_2(1, 0));
    Ray_2 r2(Point_2(0, 1), Point_2(1, 0));
    test("Ray-Ray: parallel same direction", ray_ray_intersect(r1, r2));
}

void test_ray_ray_collinear_opposite() {
    Ray_2 r1(Point_2(0, 0), Point_2(1, 0));
    Ray_2 r2(Point_2(2, 0), Point_2(-1, 0));
    test("Ray-Ray: collinear opposite", ray_ray_intersect(r1, r2));
}

void test_ray_ray_overlapping() {
    Ray_2 r1(Point_2(0, 0), Point_2(1, 0));
    Ray_2 r2(Point_2(0.5, 0), Point_2(1, 0));
    test("Ray-Ray: overlapping same direction", ray_ray_intersect(r1, r2));
}

// ============================================================================
// RAY-SEGMENT TESTS
// ============================================================================

void test_ray_segment_intersecting() {
    Ray_2 r(Point_2(0, 0), Point_2(1, 1));
    Segment_2 s(Point_2(-1, 1), Point_2(1, -1));
    test("Ray-Segment: diagonal intersect", ray_segment_intersect(r, s));
}

void test_ray_segment_parallel() {
    Ray_2 r(Point_2(0, 0), Point_2(1, 0));
    Segment_2 s(Point_2(0, 1), Point_2(1, 1));
    test("Ray-Segment: parallel", !ray_segment_intersect(r, s));
}

void test_ray_segment_endpoint() {
    Ray_2 r(Point_2(0, 0), Point_2(1, 0));
    Segment_2 s(Point_2(1, 0), Point_2(2, 1));
    test("Ray-Segment: through endpoint", ray_segment_intersect(r, s));
}

void test_ray_segment_behind() {
    Ray_2 r(Point_2(0, 0), Point_2(-1, 0));
    Segment_2 s(Point_2(1, -1), Point_2(1, 1));
    test("Ray-Segment: behind source", !ray_segment_intersect(r, s));
}

void test_ray_segment_far_away() {
    Ray_2 r(Point_2(0, 0), Point_2(1, 0));
    Segment_2 s(Point_2(10, 10), Point_2(20, 20));
    test("Ray-Segment: far away", !ray_segment_intersect(r, s));
}

void test_ray_segment_collinear() {
    Ray_2 r(Point_2(0, 0), Point_2(1, 0));
    Segment_2 s(Point_2(0.5, 0), Point_2(2, 0));
    test("Ray-Segment: collinear overlapping", ray_segment_intersect(r, s));
}

// ============================================================================
// SEGMENT-SEGMENT TESTS
// ============================================================================

void test_segment_segment_intersecting() {
    Segment_2 s1(Point_2(0, 0), Point_2(2, 2));
    Segment_2 s2(Point_2(0, 2), Point_2(2, 0));
    test("Segment-Segment: X intersection", segment_segment_intersect(s1, s2));
}

void test_segment_segment_parallel() {
    Segment_2 s1(Point_2(0, 0), Point_2(1, 0));
    Segment_2 s2(Point_2(0, 1), Point_2(1, 1));
    test("Segment-Segment: parallel", !segment_segment_intersect(s1, s2));
}

void test_segment_segment_non_intersecting() {
    Segment_2 s1(Point_2(0, 0), Point_2(1, 0));
    Segment_2 s2(Point_2(2, 0), Point_2(3, 0));
    test("Segment-Segment: separated", !segment_segment_intersect(s1, s2));
}

// ============================================================================
// MAIN TEST RUNNER
// ============================================================================

void print_summary() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "TEST SUMMARY\n";
    std::cout << std::string(70, '=') << "\n";

    int passed = 0, failed = 0;
    for (const auto& r : results) {
        if (r.passed) passed++;
        else failed++;
    }

    std::cout << "Total: " << results.size() << " tests\n";
    std::cout << "Passed: " << passed << "\n";
    std::cout << "Failed: " << failed << "\n";

    if (failed > 0) {
        std::cout << "\nFailed tests:\n";
        for (const auto& r : results) {
            if (!r.passed) {
                std::cout << "  ✗ " << r.name << "\n";
            }
        }
    }
    std::cout << std::string(70, '=') << "\n";
}

int main() {
    std::cout << "ROBUST RAY-RAY AND RAY-SEGMENT INTERSECTION TESTS\n";
    std::cout << "================================================\n\n";

    // RAY-RAY TESTS
    std::cout << "RAY-RAY INTERSECTIONS:\n";
    std::cout << std::string(70, '-') << "\n";
    test_ray_ray_same_source();
    test_ray_ray_non_intersecting();
    test_ray_ray_intersecting();
    test_ray_ray_parallel();
    test_ray_ray_collinear_opposite();
    test_ray_ray_overlapping();

    // RAY-SEGMENT TESTS
    std::cout << "\nRAY-SEGMENT INTERSECTIONS:\n";
    std::cout << std::string(70, '-') << "\n";
    test_ray_segment_intersecting();
    test_ray_segment_parallel();
    test_ray_segment_endpoint();
    test_ray_segment_behind();
    test_ray_segment_far_away();
    test_ray_segment_collinear();

    // SEGMENT-SEGMENT TESTS
    std::cout << "\nSEGMENT-SEGMENT INTERSECTIONS:\n";
    std::cout << std::string(70, '-') << "\n";
    test_segment_segment_intersecting();
    test_segment_segment_parallel();
    test_segment_segment_non_intersecting();

    // SUMMARY
    print_summary();

    // Return success if all tests passed
    int failed_count = 0;
    for (const auto& r : results) {
        if (!r.passed) failed_count++;
    }

    return (failed_count == 0) ? 0 : 1;
}
