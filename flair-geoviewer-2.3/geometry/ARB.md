An ARB is a portion of space bounded by 4, 5 or 6 plane faces.
      Its use is rather complicated and is now superseded by the
      availability of infinite planes (XYP, XZP, YZP and PLA).
      For completeness, however, a description of input will be reported
      here.
      Assign an index (1 to 8) to each vertex. For each vertex, give
      the x, y, z coordinates on 4 lines (8 lines in high-accuracy body
      fixed format), with 6 numbers on each line (3 numbers per line in
      high-accuracy fixed format):
      V1_x, V1_y, V1_z, V2_x, V2_y, V2_z,
      V3_x, V3_y, V3_z, V4_x, V4_y, V4_z,
      V5_x, V5_y, V5_z, V6_x, V6_y, V6_z,
      V7_x, V7_y, V7_z, V8_x, V8_y, V8_z
      The vertices not appearing in face descriptions are ignored, but
      eight must always be supplied. The above vertex lines are followed
      by one line describing the faces (two lines in high-accuracy fixed format).
      Each face is described by a four-digit number in floating point format,
      such that each digit gives the indices of the three or four vertex
      points of that face. These indices must be entered in the same order
      for each face (i.e. all clockwise or all counterclockwise).
      When the number of the faces is less than 6, the remaining face
      description(s) must be zero, and must appear at the end of the
      list. If a face has three vertices, the omitted position may be
      either 0 or a repeat of one of the other vertices.
      An ARB definition extends over 5 lines in default fixed format, or over
      10 lines in high-accuracy body fixed format.

      Example in default fixed format:
 *...+....1....+....2....+....3....+....4....+....5....+....6....+....7..
   ARB    5     -44.2     -36.5    3572.0     -44.2     -33.5    3572.0
                -37.1     -31.0    3572.0     -37.1     -39.0    3572.0
                -44.2     -36.5       0.0     -44.2     -33.5       0.0
                -37.1     -31.0       0.0     -37.1     -39.0       0.0
                1234.     1562.     5876.     1485.     4378.     6732.
 * (a right prism of trapezoidal base, of height 3572 cm parallel to the
 *  z axis)

      The same example in high-accuracy body fixed format:
 *...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+.
   ARB    5                 -44.2                 -36.5                3572.0
                            -44.2                 -33.5                3572.0
                            -37.1                 -31.0                3572.0
                            -37.1                 -39.0                3572.0
                            -44.2                 -36.5                   0.0
                            -44.2                 -33.5                   0.0
                            -37.1                 -31.0                   0.0
                            -37.1                 -39.0                   0.0
                             1234.                1562.                 5876.
                             1485.                4378.                 6732.

      The same, in free format:
 ARB oddbody -44.2 -36.5 3572.0 -44.2 -33.5 3572.0 -37.1 -31.0 3572.0 -37.1
          -39.0 3572.0 -44.2 -36.5 0.0 -44.2 -33.5 0.0 -37.1 -31.0 0.0 -37.1
          -39.0 0.0 1234. 1562. 5876. 1485. 4378. 6732.
