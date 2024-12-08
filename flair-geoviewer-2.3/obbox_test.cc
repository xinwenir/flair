#include "minunit.h"
#include "vector.h"
#include "obbox.h"
#include "matrix4.h"

/* Constructors */
MU_TEST(test_default) {
	OBBox bb;
	mu_assert_double_eq( INFINITE, bb.low().x);
	mu_assert_double_eq( INFINITE, bb.low().y);
	mu_assert_double_eq( INFINITE, bb.low().z);
	mu_assert_double_eq(-INFINITE, bb.high().x);
	mu_assert_double_eq(-INFINITE, bb.high().y);
	mu_assert_double_eq(-INFINITE, bb.high().z);
}

MU_TEST(test_default_not_isValid) {
	OBBox bb;
	mu_assert(!bb.isValid(), "shold not be valid");
}

MU_TEST(test_copy_constructor) {
	OBBox b;
	b.P = Vector(10.0, 20.0, 30.0);
	b.add( -1.0, -2.0, -3.0);
	b.add(  1.0,  2.0,  3.0);

	OBBox a(b);

	mu_assert_double_eq(10.0, a.P.x);
	mu_assert_double_eq(20.0, a.P.y);
	mu_assert_double_eq(30.0, a.P.z);

	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();
	mu_assert_double_eq(-1.0, low.x);
	mu_assert_double_eq(-2.0, low.y);
	mu_assert_double_eq(-3.0, low.z);
	mu_assert_double_eq( 1.0, high.x);
	mu_assert_double_eq( 2.0, high.y);
	mu_assert_double_eq( 3.0, high.z);
}

/* add */
MU_TEST(test_add_one) {
	OBBox bb;
	bb.add(-15.0, -16.0, -17.0);

	mu_assert(!bb.isValid(), "One point should not define a valid bb");
}

MU_TEST(test_add) {
	OBBox bb;
	bb.add(-15.0, -16.0, -17.0);
	bb.add(15.0, 16.0, 17.0);

	mu_assert_double_eq(-15.0, bb.low().x);
	mu_assert_double_eq(-16.0, bb.low().y);
	mu_assert_double_eq(-17.0, bb.low().z);
	mu_assert_double_eq( 15.0, bb.high().x);
	mu_assert_double_eq( 16.0, bb.high().y);
	mu_assert_double_eq( 17.0, bb.high().z);
	mu_assert(bb.isValid(), "Should be valid");
}

MU_TEST(test_add_diff_P) {
	OBBox bb;
	bb.P = Vector(10.0, 20.0, 30.0);
	bb.add(-15.0, -16.0, -17.0);
	bb.add(15.0, 16.0, 17.0);

	mu_assert_double_eq(-15.0-10.0, bb.low().x);
	mu_assert_double_eq(-16.0-20.0, bb.low().y);
	mu_assert_double_eq(-17.0-30.0, bb.low().z);
	mu_assert_double_eq( 15.0-10.0, bb.high().x);
	mu_assert_double_eq( 16.0-20.0, bb.high().y);
	mu_assert_double_eq( 17.0-30.0, bb.high().z);
	mu_assert(bb.isValid(), "Should be valid");
}

/* assignment */
MU_TEST(test_assignment) {
	OBBox a;
	a.add(-4.0, -4.0, -4.0);
	a.add( 4.0,  4.0,  4.0);

	OBBox b;
	b.P = Vector(10.0, 20.0, 30.0);
	b.add( -1.0, -2.0, -3.0);
	b.add(  1.0,  2.0,  3.0);

	a = b;

	mu_assert_double_eq(10.0, a.P.x);
	mu_assert_double_eq(20.0, a.P.y);
	mu_assert_double_eq(30.0, a.P.z);

	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();
	mu_assert_double_eq(-1.0, low.x);
	mu_assert_double_eq(-2.0, low.y);
	mu_assert_double_eq(-3.0, low.z);
	mu_assert_double_eq( 1.0, high.x);
	mu_assert_double_eq( 2.0, high.y);
	mu_assert_double_eq( 3.0, high.z);
}

/* invalidate */
MU_TEST(test_invalidate) {
	OBBox bb;
	bb.add(-15.0, -16.0, -17.0);
	bb.add(15.0, 16.0, 17.0);

	bb.invalidate();
	mu_assert(!bb.isValid(), "Should not be valid");
}


/* inside */
MU_TEST(test_inside) {
	OBBox bb;
	bb.X = Vector(0.0, 1.0, 0.0);
	bb.Y = Vector(-1.0, 0.0, 0.0);
	bb.Z = Vector(0.0, 0.0, 1.0);
	bb.add(-15.0, -1.0,-1.0);
	bb.add( 15.0,  1.0, 1.0);

	mu_assert(bb.inside( 0.0,  0.0, 0.0), "Should be inside");
	mu_assert(bb.inside(-14.5, 0.9, 0.9), "Should be inside");
	mu_assert(bb.inside( 14.5, 0.9, 0.9), "Should be inside");

	mu_assert(!bb.inside(-15.5, 1.5, 0.0), "Should be outside");
	mu_assert(!bb.inside( 15.5, 1.5, 0.0), "Should be outside");
	mu_assert(!bb.inside(  0.0, 1.5, 0.0), "Should be outside");
}

MU_TEST(test_inside_diff_P) {
	OBBox bb;
	bb.X = Vector(0.0, 1.0, 0.0);
	bb.Y = Vector(-1.0, 0.0, 0.0);
	bb.Z = Vector(0.0, 0.0, 1.0);
	bb.P = Vector(10.0, 20.0, 30.0);
	bb.add(-15.0, -1.0,-1.0);
	bb.add( 15.0,  1.0, 1.0);

	mu_assert(bb.inside( 0.0,  0.0, 0.0), "Should be inside");
	mu_assert(bb.inside(-14.5, 0.9, 0.9), "Should be inside");
	mu_assert(bb.inside( 14.5, 0.9, 0.9), "Should be inside");

	mu_assert(!bb.inside(-15.5, 1.5, 0.0), "Should be outside");
	mu_assert(!bb.inside( 15.5, 1.5, 0.0), "Should be outside");
	mu_assert(!bb.inside(  0.0, 1.5, 0.0), "Should be outside");
}

MU_TEST(test_inside_boundary) {
	OBBox bb;
	bb.X = Vector(0.0, 1.0, 0.0);
	bb.Y = Vector(-1.0, 0.0, 0.0);
	bb.Z = Vector(0.0, 0.0, 1.0);
	bb.add(-15.0, -1.0,-1.0);
	bb.add( 15.0,  1.0, 1.0);

	mu_assert(bb.inside( 0.0,  0.0, 0.0), "Should be inside");
	mu_assert(bb.inside(-14.5, 0.9, 0.9), "Should be inside");
	mu_assert(bb.inside( 14.5, 0.9, 0.9), "Should be inside");

	mu_assert(!bb.inside(-15.5, 1.5, 0.0), "Should be outside");
	mu_assert(!bb.inside( 15.5, 1.5, 0.0), "Should be outside");
	mu_assert(!bb.inside(  0.0, 1.5, 0.0), "Should be outside");
}

MU_TEST(test_inside_boundary_diff_P) {
	OBBox bb;
	bb.X = Vector(0.0, 1.0, 0.0);
	bb.Y = Vector(-1.0, 0.0, 0.0);
	bb.Z = Vector(0.0, 0.0, 1.0);
	bb.P = Vector(10.0, 20.0, 30.0);
	bb.add(-15.0, -1.0,-1.0);
	bb.add( 15.0,  1.0, 1.0);

	mu_assert(bb.inside( 0.0,  0.0, 0.0), "Should be inside");
	mu_assert(bb.inside(-14.5, 0.9, 0.9), "Should be inside");
	mu_assert(bb.inside( 14.5, 0.9, 0.9), "Should be inside");

	mu_assert(!bb.inside(-15.5, 1.5, 0.0), "Should be outside");
	mu_assert(!bb.inside( 15.5, 1.5, 0.0), "Should be outside");
	mu_assert(!bb.inside(  0.0, 1.5, 0.0), "Should be outside");
}

/* vertex */
MU_TEST(test_vertex) {
	OBBox bb;
	bb.add(-15.0, -16.0, -17.0);
	bb.add( 15.0,  16.0,  17.0);

	Vector low  =  bb.lowVertex();
	Vector high = bb.highVertex();

	mu_assert_double_eq(-15.0, low.x);
	mu_assert_double_eq(-16.0, low.y);
	mu_assert_double_eq(-17.0, low.z);

	mu_assert_double_eq( 15.0, high.x);
	mu_assert_double_eq( 16.0, high.y);
	mu_assert_double_eq( 17.0, high.z);
}

MU_TEST(test_vertex_diff_P) {
	OBBox bb;
	bb.P = Vector(10.0, 20.0, 30.0);
	bb.add(-15.0, -16.0, -17.0);
	bb.add( 15.0,  16.0,  17.0);

	Vector low  =  bb.lowVertex();
	Vector high = bb.highVertex();

	mu_assert_double_eq(-15.0, low.x);
	mu_assert_double_eq(-16.0, low.y);
	mu_assert_double_eq(-17.0, low.z);

	mu_assert_double_eq( 15.0, high.x);
	mu_assert_double_eq( 16.0, high.y);
	mu_assert_double_eq( 17.0, high.z);
}

/* Intersect */
MU_TEST(test_intersect_SameCoordSystem_equal) {
	OBBox a;
	a.add(-4.0, -4.0, -4.0);
	a.add( 4.0,  4.0,  4.0);

	OBBox b;
	b.add(-4.0, -4.0, -4.0);
	b.add( 4.0,  4.0,  4.0);

	a.Intersect(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert_double_eq(-4.0, low.x);
	mu_assert_double_eq(-4.0, low.y);
	mu_assert_double_eq(-4.0, low.z);
	mu_assert_double_eq( 4.0, high.x);
	mu_assert_double_eq( 4.0, high.y);
	mu_assert_double_eq( 4.0, high.z);
}

MU_TEST(test_intersect_SameCoordSystem_inside) {
	OBBox a;
	a.add(-4.0, -4.0, -4.0);
	a.add( 4.0,  4.0,  4.0);

	OBBox b;
	b.add(-2.0, -2.0, -4.0);
	b.add( 2.0,  2.0,  4.0);

	a.Intersect(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert_double_eq(-2.0, low.x);
	mu_assert_double_eq(-2.0, low.y);
	mu_assert_double_eq(-4.0, low.z);
	mu_assert_double_eq( 2.0, high.x);
	mu_assert_double_eq( 2.0, high.y);
	mu_assert_double_eq( 4.0, high.z);
}

MU_TEST(test_intersect_SameCoordSystem_through) {
	OBBox a;
	a.add(-4.0, -4.0, -4.0);
	a.add( 4.0,  4.0,  4.0);

	OBBox b;
	b.add(-4.0, -4.0, -6.0);
	b.add( 2.0,  2.0,  6.0);

	a.Intersect(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert_double_eq(-4.0, low.x);
	mu_assert_double_eq(-4.0, low.y);
	mu_assert_double_eq(-4.0, low.z);
	mu_assert_double_eq( 2.0, high.x);
	mu_assert_double_eq( 2.0, high.y);
	mu_assert_double_eq( 4.0, high.z);
}


/* Same coordinate system different origin */
MU_TEST(test_intersect_SameCoordSystemDiffP_equal1) {
	OBBox a;
	a.add(-4.0, -4.0, -4.0);
	a.add( 4.0,  4.0,  4.0);

	OBBox b;
	b.P = Vector(10.0, 20.0, 30.0);
	b.add(-4.0, -4.0, -4.0);
	b.add(  4.0,  4.0,  4.0);

	a.Intersect(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert(a.isValid(), "Result should be valid");
	mu_assert_double_eq(-4.0, low.x);
	mu_assert_double_eq(-4.0, low.y);
	mu_assert_double_eq(-4.0, low.z);
	mu_assert_double_eq( 4.0, high.x);
	mu_assert_double_eq( 4.0, high.y);
	mu_assert_double_eq( 4.0, high.z);
}
MU_TEST(test_intersect_SameCoordSystemDiffP_equal2) {
	OBBox a;
	a.P = Vector(10.0, 20.0, 30.0);
	a.add(-4.0, -4.0, -4.0);
	a.add( 4.0,  4.0,  4.0);

	OBBox b;
	b.add(-4.0, -4.0, -4.0);
	b.add(  4.0,  4.0,  4.0);

	a.Intersect(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert(a.isValid(), "Result should be valid");
	mu_assert_double_eq(-4.0, low.x);
	mu_assert_double_eq(-4.0, low.y);
	mu_assert_double_eq(-4.0, low.z);
	mu_assert_double_eq( 4.0, high.x);
	mu_assert_double_eq( 4.0, high.y);
	mu_assert_double_eq( 4.0, high.z);
}

MU_TEST(test_intersect_SameCoordSystemDiffP_inside) {
	OBBox a;
	a.add(-4.0, -4.0, -4.0);
	a.add( 4.0,  4.0,  4.0);

	OBBox b;
	Vector Pb(10.0, 20.0, 30.0);
	b.P = Pb;
	b.add(-2.0, -2.0, -4.0);
	b.add( 2.0,  2.0,  4.0);

	a.Intersect(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert_double_eq(-2.0, low.x);
	mu_assert_double_eq(-2.0, low.y);
	mu_assert_double_eq(-4.0, low.z);
	mu_assert_double_eq( 2.0, high.x);
	mu_assert_double_eq( 2.0, high.y);
	mu_assert_double_eq( 4.0, high.z);
}

MU_TEST(test_intersect_SameCoordSystemDiffP_through) {
	OBBox a;
	a.add(-4.0, -4.0, -4.0);
	a.add( 4.0,  4.0,  4.0);

	OBBox b;
	Vector Pb(10.0, 20.0, 30.0);
	b.P = Pb;
	b.add(-2.0, -2.0, -6.0);
	b.add( 2.0,  2.0, 6.0);

	a.Intersect(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert_double_eq(-2.0, low.x);
	mu_assert_double_eq(-2.0, low.y);
	mu_assert_double_eq(-4.0, low.z);
	mu_assert_double_eq( 2.0, high.x);
	mu_assert_double_eq( 2.0, high.y);
	mu_assert_double_eq( 4.0, high.z);
}

/* Different coordinate system + different origin */
MU_TEST(test_intersect_DiffCoordSystem_equal) {
	OBBox a;
	a.add(-4.0, -5.0, -6.0);
	a.add( 4.0,  5.0,  6.0);

	OBBox b;
	Vector Pb(10.0, 20.0, 30.0);
	b.P = Pb;
	b.X = Vector::Yo;
	b.Y = Vector::Xo;
	b.Z = Vector::Zo;
	b.add(-4.0, -5.0, -6.0);
	b.add( 4.0,  5.0,  6.0);

	a.Intersect(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert_double_eq(-4.0, low.x);
	mu_assert_double_eq(-5.0, low.y);
	mu_assert_double_eq(-6.0, low.z);
	mu_assert_double_eq( 4.0, high.x);
	mu_assert_double_eq( 5.0, high.y);
	mu_assert_double_eq( 6.0, high.z);
}

MU_TEST(test_intersect_DiffCoordSystem_inside) {
	OBBox a;
	a.add(-4.0, -5.0, -6.0);
	a.add( 4.0,  5.0,  6.0);

	OBBox b;
	Vector Pb(10.0, 20.0, 30.0);
	b.P = Pb;
	b.X = Vector::Yo;
	b.Y = Vector::Xo;
	b.Z = Vector::Zo;
	b.add(-2.0, -3.0, -4.0);
	b.add( 2.0,  3.0,  4.0);

	a.Intersect(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert_double_eq(-2.0, low.x);
	mu_assert_double_eq(-3.0, low.y);
	mu_assert_double_eq(-4.0, low.z);
	mu_assert_double_eq( 2.0, high.x);
	mu_assert_double_eq( 3.0, high.y);
	mu_assert_double_eq( 4.0, high.z);
}

MU_TEST(test_intersect_DiffCoordSystem_through) {
	OBBox a;
	a.add(-4.0, -5.0, -6.0);
	a.add( 4.0,  5.0,  6.0);

	OBBox b;
	Vector Pb(10.0, 20.0, 30.0);
	b.P = Pb;
	b.X = Vector::Yo;
	b.Y = Vector::Xo;
	b.Z = Vector::Zo;
	b.add(-2.0, -3.0, -8.0);
	b.add( 2.0,  3.0,  8.0);

	a.Intersect(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert_double_eq(-2.0, low.x);
	mu_assert_double_eq(-3.0, low.y);
	mu_assert_double_eq(-6.0, low.z);
	mu_assert_double_eq( 2.0, high.x);
	mu_assert_double_eq( 3.0, high.y);
	mu_assert_double_eq( 6.0, high.z);
}

MU_TEST(test_intersect_DiffCoordSystem_insideRotated) {
	OBBox a;
	a.add(-4.0, -4.0, -4.0);
	a.add( 4.0,  4.0,  4.0);

	OBBox b;
	b.X = Vector(0.5, 0.5, 0.0);
	b.X.normalize();
	b.Y = Vector(-0.5, 0.5, 0.0);
	b.Y.normalize();
	b.Z = Vector::Zo;
	b.low() = Vector(-2.0, -2.0, -2.0);
	b.high() = Vector(2.0,  2.0,  2.0);

	a.Intersect(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert_double_eq(             0.0, low.x);
	mu_assert_double_eq(-2.0 * sqrt(2.0), low.y);
	mu_assert_double_eq(            -2.0, low.z);

	mu_assert_double_eq(             0.0, high.x);
	mu_assert_double_eq( 2.0 * sqrt(2.0), high.y);
	mu_assert_double_eq(             2.0, high.z);
}

/* Intersect */
MU_TEST(test_intersect_ReturnsSmallestVolume_bSmaller) {
	OBBox a;
	a.add(-4.0, -4.0, -4.0);
	a.add( 4.0,  4.0,  4.0);

	OBBox b;
	b.X = Vector(0.5, 0.5, 0.0);
	b.X.normalize();
	b.Y = Vector(-0.5, 0.5, 0.0);
	b.Y.normalize();
	b.Z = Vector::Zo;
	b.low() = Vector(-2.0, -2.0, -2.0);
	b.high() = Vector( 2.0,  2.0,  2.0);

	a.Intersect(b);

	mu_assert_double_eq(64.0, a.volume());
}

MU_TEST(test_union) {
	OBBox a;
	a.add(-4.0, -4.0, -4.0);
	a.add( 4.0,  4.0,  4.0);

	OBBox b;
	b.add(-4.0, -4.0, 4.0);
	b.add( 4.0,  4.0, 8.0);

	a.Union(b);

	mu_assert_double_eq(128.0, a.volume());
}

/* Infinite bodies */
MU_TEST(test_intersect_infBodies1) {
	OBBox a;
	// Infinite along x
	a.low().set(-INFINITE, -2.0, -2.0);
	a.high().set( INFINITE,  2.0,  2.0);

	OBBox b;
	// Infinite along 3 dimensions but with a lower limit in x
	b.low().set( 0.0, -INFINITE, -INFINITE);
	b.high().set( INFINITE,  INFINITE,  INFINITE);

	a.Intersect(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert_double_eq( 0.0, low.x);
	mu_assert_double_eq(-2.0, low.y);
	mu_assert_double_eq(-2.0, low.z);
	mu_assert_double_eq( INFINITE, high.x);
	mu_assert_double_eq( 2.0, high.y);
	mu_assert_double_eq( 2.0, high.z);
}

MU_TEST(test_intersect_infBodies2) {
	OBBox a;
	// Infinite along x
	a.low().set(-INFINITE, -2.0, -2.0);
	a.high().set( INFINITE,  2.0,  2.0);

	OBBox b;
	// Infinite along 3 dimensions but with a higher limit in x
	b.low().set( -INFINITE, -INFINITE, -INFINITE);
	b.high().set( 0.0,  INFINITE,  INFINITE);

	a.Intersect(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert_double_eq(-INFINITE, low.x);
	mu_assert_double_eq(-2.0, low.y);
	mu_assert_double_eq(-2.0, low.z);
	mu_assert_double_eq( 0.0, high.x);
	mu_assert_double_eq( 2.0, high.y);
	mu_assert_double_eq( 2.0, high.z);
}

MU_TEST(test_intersect_infBodies3) {
	OBBox a;
	// Infinite along x
	a.low().set(-2.0, -2.0, -2.0);
	a.high().set(2.0,  2.0,  2.0);

	OBBox b;
	// Infinite along 3 dimensions but with a higher limit in x
	b.low().set( -INFINITE, 0.0, -INFINITE);
	b.high().set( INFINITE,  INFINITE,  INFINITE);

	a.Intersect(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert_double_eq(-2.0, low.x);
	mu_assert_double_eq( 0.0, low.y);
	mu_assert_double_eq(-2.0, low.z);
	mu_assert_double_eq( 2.0, high.x);
	mu_assert_double_eq( 2.0, high.y);
	mu_assert_double_eq( 2.0, high.z);
}

MU_TEST(test_difference1) {
	OBBox a;
	a.low().set(  0.0,-5.0, 0.0);
	a.high().set( 5.0, 0.0, 5.0);

	OBBox b;
	b.X.set( 0.5,  0.5, 0.0);
	b.Y.set( 0.5, -0.5, 0.0);
	b.Z.set( 0.0,  0.0, 0.0);
	b.P.set( 2.0, -5.0, 0.0);
	b.low().set( -2.0, -4.0, -1.0);
	b.high().set( 6.0,  0.0,  6.0);

	a.Difference(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert_double_eq( 0.0, low.x);
	mu_assert_double_eq(-5.0, low.y);
	mu_assert_double_eq( 0.0, low.z);
	mu_assert_double_eq( 5.0, high.x);
	mu_assert_double_eq( 0.0, high.y);
	mu_assert_double_eq( 5.0, high.z);
}

MU_TEST(test_difference_infBodies1) {
	OBBox a;
	// Infinite along x
	a.low().set(-INFINITE, -2.0, -2.0);
	a.high().set( INFINITE,  2.0,  2.0);

	OBBox b;
	// Infinite along 3 dimensions but with a lower limit in x
	b.low().set( 0.0, -INFINITE, -INFINITE);
	b.high().set( INFINITE,  INFINITE,  INFINITE);

	a.Difference(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert_double_eq(-INFINITE, low.x);
	mu_assert_double_eq(-2.0, low.y);
	mu_assert_double_eq(-2.0, low.z);
	mu_assert_double_eq( 0.0, high.x);
	mu_assert_double_eq( 2.0, high.y);
	mu_assert_double_eq( 2.0, high.z);
}

MU_TEST(test_difference_infBodies2) {
	OBBox a;
	// Infinite along x
	a.low().set(-INFINITE, -2.0, -2.0);
	a.high().set( INFINITE,  2.0,  2.0);

	OBBox b;
	// Infinite along 3 dimensions but with a lower limit in x
	b.low().set( -INFINITE, -INFINITE, -INFINITE);
	b.high().set( 0.0,  INFINITE,  INFINITE);

	a.Difference(b);
	Vector low  =  a.lowVertex();
	Vector high = a.highVertex();

	mu_assert_double_eq( 0.0, low.x);
	mu_assert_double_eq(-2.0, low.y);
	mu_assert_double_eq(-2.0, low.z);
	mu_assert_double_eq( INFINITE, high.x);
	mu_assert_double_eq( 2.0, high.y);
	mu_assert_double_eq( 2.0, high.z);
}

MU_TEST(test_transform_unit) {
	OBBox bb;
	bb.add(-15.0, -16.0, -17.0);
	bb.add(15.0, 16.0, 17.0);

	Matrix4 m;
	m(0,0) = 1.0; m(0,1) = 0.0; m(0,2) = 0.0; m(0,3) = 0.0;
	m(1,0) = 0.0; m(1,1) = 1.0; m(1,2) = 0.0; m(1,3) = 0.0;
	m(2,0) = 0.0; m(2,1) = 0.0; m(2,2) = 1.0; m(2,3) = 0.0;
	m(3,0) = 0.0; m(3,1) = 0.0; m(3,2) = 0.0; m(3,3) = 1.0;

	mu_assert_double_eq(32640.0, bb.volume());
	bb.transform(m);
	mu_assert_double_eq(-15.0, bb.low().x);
	mu_assert_double_eq(-16.0, bb.low().y);
	mu_assert_double_eq(-17.0, bb.low().z);
	mu_assert_double_eq( 15.0, bb.high().x);
	mu_assert_double_eq( 16.0, bb.high().y);
	mu_assert_double_eq( 17.0, bb.high().z);
	mu_assert(bb.isValid(), "Should be valid");
	mu_assert_double_eq(32640.0, bb.volume());
}

MU_TEST(test_transform_translate) {
	OBBox bb;
	bb.add(-15.0, -16.0, -17.0);
	bb.add(15.0, 16.0, 17.0);

	Matrix4 m;
	m(0,0) = 1.0; m(0,1) = 0.0; m(0,2) = 0.0; m(0,3) = 10.0;
	m(1,0) = 0.0; m(1,1) = 1.0; m(1,2) = 0.0; m(1,3) = 20.0;
	m(2,0) = 0.0; m(2,1) = 0.0; m(2,2) = 1.0; m(2,3) = 30.0;
	m(3,0) = 0.0; m(3,1) = 0.0; m(3,2) = 0.0; m(3,3) = 1.0;

	mu_assert_double_eq(32640.0, bb.volume());
	bb.transform(m);
	mu_assert_double_eq(-15.0, bb.lowVertex().x);
	mu_assert_double_eq(-16.0, bb.lowVertex().y);
	mu_assert_double_eq(-17.0, bb.lowVertex().z);
	mu_assert_double_eq( 15.0, bb.highVertex().x);
	mu_assert_double_eq( 16.0, bb.highVertex().y);
	mu_assert_double_eq( 17.0, bb.highVertex().z);
	mu_assert(bb.isValid(), "Should be valid");
	mu_assert_double_eq(32640.0, bb.volume());
}

MU_TEST(test_transform_rotate) {
	OBBox bb;
	bb.add(-15.0, -16.0, -17.0);
	bb.add(15.0, 16.0, 17.0);

	Matrix4 m;
	m(0,0) = 0.0; m(0,1) =-1.0; m(0,2) = 0.0; m(0,3) = 0.0;
	m(1,0) = 1.0; m(1,1) = 0.0; m(1,2) = 0.0; m(1,3) = 0.0;
	m(2,0) = 0.0; m(2,1) = 0.0; m(2,2) = 1.0; m(2,3) = 0.0;
	m(3,0) = 0.0; m(3,1) = 0.0; m(3,2) = 0.0; m(3,3) = 1.0;

	mu_assert_double_eq(32640.0, bb.volume());
	bb.transform(m);
	mu_assert_double_eq(15.0, Abs(bb.lowVertex().x));
	mu_assert_double_eq(16.0, Abs(bb.lowVertex().y));
	mu_assert_double_eq(17.0, Abs(bb.lowVertex().z));
	mu_assert_double_eq(15.0, Abs(bb.highVertex().x));
	mu_assert_double_eq(16.0, Abs(bb.highVertex().y));
	mu_assert_double_eq(17.0, Abs(bb.highVertex().z));
	mu_assert(bb.isValid(), "Should be valid");
	mu_assert_double_eq(32640.0, bb.volume());
}

MU_TEST(test_transform_unit_P) {
	OBBox bb;
	bb.P.set(10.0, 20.0, 30.0);
	bb.add(-15.0, -16.0, -17.0);
	bb.add(15.0, 16.0, 17.0);

	Matrix4 m;
	m(0,0) = 1.0; m(0,1) = 0.0; m(0,2) = 0.0; m(0,3) = 0.0;
	m(1,0) = 0.0; m(1,1) = 1.0; m(1,2) = 0.0; m(1,3) = 0.0;
	m(2,0) = 0.0; m(2,1) = 0.0; m(2,2) = 1.0; m(2,3) = 0.0;
	m(3,0) = 0.0; m(3,1) = 0.0; m(3,2) = 0.0; m(3,3) = 1.0;

	Vector low  = bb.lowVertex();
	Vector high = bb.highVertex();

	mu_assert_double_eq(32640.0, bb.volume());
	bb.transform(m);
	mu_assert_double_eq(-15.0, low.x);
	mu_assert_double_eq(-16.0, low.y);
	mu_assert_double_eq(-17.0, low.z);
	mu_assert_double_eq( 15.0, high.x);
	mu_assert_double_eq( 16.0, high.y);
	mu_assert_double_eq( 17.0, high.z);
	mu_assert(bb.isValid(), "Should be valid");
	mu_assert_double_eq(32640.0, bb.volume());
}

MU_TEST(test_transform_translate_P) {
	OBBox bb;
	bb.P.set(10.0, 20.0, 30.0);
	bb.add(-15.0, -16.0, -17.0);
	bb.add(15.0, 16.0, 17.0);

	Matrix4 m;
	m(0,0) = 1.0; m(0,1) = 0.0; m(0,2) = 0.0; m(0,3) = 10.0;
	m(1,0) = 0.0; m(1,1) = 1.0; m(1,2) = 0.0; m(1,3) = 20.0;
	m(2,0) = 0.0; m(2,1) = 0.0; m(2,2) = 1.0; m(2,3) = 30.0;
	m(3,0) = 0.0; m(3,1) = 0.0; m(3,2) = 0.0; m(3,3) = 1.0;

	Vector low  = bb.lowVertex();
	Vector high = bb.highVertex();

	mu_assert_double_eq(32640.0, bb.volume());
	bb.transform(m);
	mu_assert_double_eq(-15.0, low.x);
	mu_assert_double_eq(-16.0, low.y);
	mu_assert_double_eq(-17.0, low.z);
	mu_assert_double_eq( 15.0, high.x);
	mu_assert_double_eq( 16.0, high.y);
	mu_assert_double_eq( 17.0, high.z);
	mu_assert(bb.isValid(), "Should be valid");
	mu_assert_double_eq(32640.0, bb.volume());
}

MU_TEST(test_transform_rotate_P) {
	OBBox bb;
	bb.P.set(10.0, 20.0, 30.0);
	bb.add(-15.0, -16.0, -17.0);
	bb.add(15.0, 16.0, 17.0);

	Matrix4 m;
	m(0,0) = 0.0; m(0,1) =-1.0; m(0,2) = 0.0; m(0,3) = 0.0;
	m(1,0) = 1.0; m(1,1) = 0.0; m(1,2) = 0.0; m(1,3) = 0.0;
	m(2,0) = 0.0; m(2,1) = 0.0; m(2,2) = 1.0; m(2,3) = 0.0;
	m(3,0) = 0.0; m(3,1) = 0.0; m(3,2) = 0.0; m(3,3) = 1.0;

	Vector low  = bb.lowVertex();
	Vector high = bb.highVertex();

	mu_assert_double_eq(32640.0, bb.volume());
	bb.transform(m);
	mu_assert_double_eq(-15.0, low.x);
	mu_assert_double_eq(-16.0, low.y);
	mu_assert_double_eq(-17.0, low.z);
	mu_assert_double_eq( 15.0, high.x);
	mu_assert_double_eq( 16.0, high.y);
	mu_assert_double_eq( 17.0, high.z);
	mu_assert(bb.isValid(), "Should be valid");
	mu_assert_double_eq(32640.0, bb.volume());
}

MU_TEST(test_plane_intersect1) {
	OBBox bb;
	bb.P.set(10.0, 20.0, 30.0);
	bb.X = Vector(0.0, 1.0, 0.0);
	bb.Y = Vector(-1.0, 0.0, 0.0);
	bb.Z = Vector(0.0, 0.0, 1.0);
	bb.add(-15.0, -16.0, -17.0);
	bb.add(15.0, 16.0, 17.0);

	// Matrix defining the plane
	Matrix4 matrix;
	matrix.identity();
	matrix(0,3) = 5.0; matrix(1,3) = 5.0; matrix(2,3) = 5.0; // position

	// Intersect
	int location = bb.intersectWithPlane(matrix, 0.0, 50.0, 0.0, 50.0);
	// Check results
	mu_assert_int_eq(2, location);
}

MU_TEST(test_plane_intersect2) {
	OBBox bb;
	bb.P.set(10.0, 20.0, 30.0);
	bb.X = Vector(0.0, 1.0, 0.0);
	bb.Y = Vector(-1.0, 0.0, 0.0);
	bb.Z = Vector(0.0, 0.0, 1.0);
	bb.add(-15.0, -16.0, -17.0);
	bb.add(15.0, 16.0, 17.0);

	// Create plane
	Matrix4 matrix;
	matrix.identity();
	matrix(0,3) = 5.0; matrix(1,3) = 5.0; matrix(2,3) = 5.0; // position

	// Intersect
	int location = bb.intersectWithPlane(matrix, 0.0, 5.0, 0.0, 5.0);
	// Check results
	mu_assert_int_eq(2, location); // This is a bounding box, could be filled or not
}

MU_TEST(test_plane_intersect3) {
	OBBox bb;
	bb.P.set(10.0, 20.0, 30.0);
	bb.X = Vector(0.0, 1.0, 0.0);
	bb.Y = Vector(-1.0, 0.0, 0.0);
	bb.Z = Vector(0.0, 0.0, 1.0);
	bb.add(-15.0, -16.0, -17.0);
	bb.add(15.0, 16.0, 17.0);

	// Matrix defining the plane
	Matrix4 matrix;
	matrix.identity();
	matrix(0,3) = 25.0; matrix(1,3) = 25.0; matrix(2,3) = 25.0; // position

	// Intersect
	int location = bb.intersectWithPlane(matrix, 0.0, 5.0, 0.0, 5.0);
	// Check results
	mu_assert_int_eq(0, location);
}

/** Test Suite **/
MU_TEST_SUITE(test_suite) {
	MU_RUN_TEST(test_default);
	MU_RUN_TEST(test_default_not_isValid);
	MU_RUN_TEST(test_copy_constructor);

	MU_RUN_TEST(test_add_one);
	MU_RUN_TEST(test_add);
	MU_RUN_TEST(test_add_diff_P);

	MU_RUN_TEST(test_assignment);

	MU_RUN_TEST(test_invalidate);

	MU_RUN_TEST(test_vertex);
	MU_RUN_TEST(test_vertex_diff_P);

	MU_RUN_TEST(test_inside);
	MU_RUN_TEST(test_inside_diff_P);
	MU_RUN_TEST(test_inside_boundary);
	MU_RUN_TEST(test_inside_boundary_diff_P);

	MU_RUN_TEST(test_intersect_SameCoordSystem_equal);
	MU_RUN_TEST(test_intersect_SameCoordSystem_inside);
	MU_RUN_TEST(test_intersect_SameCoordSystem_through);

	MU_RUN_TEST(test_intersect_SameCoordSystemDiffP_equal1);
	MU_RUN_TEST(test_intersect_SameCoordSystemDiffP_equal2);
	MU_RUN_TEST(test_intersect_SameCoordSystemDiffP_inside);
	MU_RUN_TEST(test_intersect_SameCoordSystemDiffP_through);

	MU_RUN_TEST(test_intersect_DiffCoordSystem_equal);
	MU_RUN_TEST(test_intersect_DiffCoordSystem_inside);
	MU_RUN_TEST(test_intersect_DiffCoordSystem_through);

	MU_RUN_TEST(test_intersect_DiffCoordSystem_insideRotated);

	MU_RUN_TEST(test_intersect_ReturnsSmallestVolume_bSmaller);

	MU_RUN_TEST(test_union);

	MU_RUN_TEST(test_intersect_infBodies1);
	MU_RUN_TEST(test_intersect_infBodies2);
	MU_RUN_TEST(test_intersect_infBodies3);

	MU_RUN_TEST(test_difference1);

	MU_RUN_TEST(test_difference_infBodies1);
	MU_RUN_TEST(test_difference_infBodies2);

	MU_RUN_TEST(test_transform_unit);
	MU_RUN_TEST(test_transform_translate);
	MU_RUN_TEST(test_transform_rotate);

	MU_RUN_TEST(test_transform_unit_P);
	MU_RUN_TEST(test_transform_translate_P);
	MU_RUN_TEST(test_transform_rotate_P);

	MU_RUN_TEST(test_plane_intersect1);
	MU_RUN_TEST(test_plane_intersect2);
	MU_RUN_TEST(test_plane_intersect3);
}

int main(int, char **) {
	MU_RUN_SUITE(test_suite);

	MU_REPORT();
	return 0;
}
