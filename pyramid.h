/*****************************************************************************/
/*                                                                           */
/*  (pyramid.h)                                                              */
/*                                                                           */
/*  Include file for programs that call Pyramid.                             */
/*                                                                           */
/*  Accompanies Pyramid Version 0.46                                         */
/*  June 4, 1998                                                             */
/*                                                                           */
/*  Copyright 1998                                                           */
/*  Jonathan Richard Shewchuk                                                */
/*  Computer Sciences Division                                               */
/*  University of California at Berkeley                                     */
/*  Berkeley, California  94720-1776                                         */
/*  jrs@cs.berkeley.edu                                                      */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  How to call Pyramid from another program                                 */
/*                                                                           */
/*                                                                           */
/*  If you haven't read Pyramid's instructions (run "pyramid -h" to read     */
/*  them), you won't understand what follows.                                */
/*                                                                           */
/*  Pyramid must be compiled into an object file (pyramid.o) with the        */
/*  PYRLIBRARY symbol defined (preferably by using the -DPYRLIBRARY compiler */
/*  switch).  The makefile included with Pyramid will do this for you if     */
/*  you run "make pyrlibrary".  The resulting object file can be called via  */
/*  the procedure tetrahedralize() (if you can spell it!).                   */
/*                                                                           */
/*  If the size of the object file is important to you, you may wish to      */
/*  generate a reduced version of pyramid.o.  The CDT_ONLY symbol gets rid   */
/*  of all meshing algorithms above and beyond constrained Delaunay          */
/*  tetrahedralization.  Specifically, the -DCDT_ONLY switch eliminates      */
/*  Pyramid's -r, -q, -a, -u, and -S switches.                               */
/*                                                                           */
/*  IMPORTANT:  These definitions (PYRLIBRARY, CDT_ONLY) must be made in the */
/*  makefile or in pyramid.c itself.  Putting these definitions in this file */
/*  will not create the desired effect.                                      */
/*                                                                           */
/*                                                                           */
/*  The calling convention for tetrahedralize() follows.                     */
/*                                                                           */
/*      void tetrahedralize(pyrswitches, in, out, vorout)                    */
/*      char *pyrswitches;                                                   */
/*      struct tetrahedralizeio *in;                                         */
/*      struct tetrahedralizeio *out;                                        */
/*      struct tetrahedralizeio *vorout;                                     */
/*                                                                           */
/*  `pyrswitches' is a string containing the command line switches you wish  */
/*  to invoke.  No initial dash is required.  Some suggestions:              */
/*                                                                           */
/*  - You'll probably find it convenient to use the `z' switch so that       */
/*    points (and other items) are numbered from zero.  This simplifies      */
/*    indexing, because the first item of any type always starts at index    */
/*    [0] of the corresponding array, whether that item's number is zero or  */
/*    one.                                                                   */
/*  - You'll probably want to use the `Q' (quiet) switch in your final code, */
/*    but you can take advantage of Pyramid's printed output (including the  */
/*    `V' switch) while debugging.                                           */
/*  - If you are not using the `p', `q', `a', or `u' switches, then the      */
/*    output points will be identical to the input points, except possibly   */
/*    for the boundary markers.  If you don't need the boundary markers, you */
/*    should use the `N' (no nodes output) switch to save memory.  (If you   */
/*    do need boundary markers, but need to save memory, a good nasty trick  */
/*    is to set out->pointlist equal to in->pointlist before calling         */
/*    tetrahedralize(), so that Pyramid overwrites the input points with     */
/*    identical copies.)                                                     */
/*  - The `I' (no iteration numbers) and `g' (.off file output) switches     */
/*    have no effect when Pyramid is compiled with PYRLIBRARY defined.       */
/*                                                                           */
/*  `in', `out', and `vorout' are descriptions of the input, the output,     */
/*  and the Voronoi output.  If the `v' (Voronoi output) switch is not used, */
/*  `vorout' may be NULL.  `in' and `out' may never be NULL.                 */
/*                                                                           */
/*  Certain fields of the input and output structures must be initialized,   */
/*  as described below.                                                      */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  The `tetrahedralizeio' structure.                                        */
/*                                                                           */
/*  Used to pass data into and out of the tetrahedralize() procedure.        */
/*                                                                           */
/*                                                                           */
/*  Arrays are used to store points, tetrahedra, markers, and so forth.  In  */
/*  all cases, the first item in any array is stored starting at index [0].  */
/*  However, that item is item number `1' unless the `z' switch is used, in  */
/*  which case it is item number `0'.  Hence, you may find it easier to      */
/*  index points and segments (and tetrahedra in the neighbor list) if you   */
/*  use the `z' switch.  Unless, of course, you're calling Pyramid from a    */
/*  Fortran program.                                                         */
/*                                                                           */
/*  Description of fields (except the `numberof' fields, which are obvious): */
/*                                                                           */
/*  `pointlist':  An array of point coordinates.  The first point's x        */
/*    coordinate is at index [0], its y coordinate is at index [1], and its  */
/*    z coordinate is at index [2], followed by the coordinates of the       */
/*    remaining points.  Each point occupies three REALs.                    */
/*  `pointattributelist':  An array of point attributes.  Each point's       */
/*    attributes occupy `numberofpointattributes' REALs.                     */
/*  `pointmarkerlist':  An array of point markers; one int per point.        */
/*                                                                           */
/*  `tetrahedronlist':  An array of tetrahedron corners.  The first          */
/*    tetrahedron's first corner is at index [0], followed by its other      */
/*    three corners in the correct order, followed by any other nodes if the */
/*    tetrahedron represents a nonlinear element.  Each tetrahedron occupies */
/*    `numberofcorners' ints.                                                */
/*  `tetrahedronattributelist':  An array of tetrahedron attributes.  Each   */
/*    tetrahedron's attributes occupy `numberoftetrahedronattributes' REALs. */
/*  `tetrahedronvolumelist':  An array of tetrahedron volume constraints;    */
/*    one REAL per tetrahedron.  Input only.                                 */
/*  `neighborlist':  An array of tetrahedron neighbors; four ints per        */
/*    tetrahedron.  Output only.                                             */
/*                                                                           */
/*  `segmentlist':  An array of segment endpoints.  The first segment's      */
/*    endpoints are at indices [0] and [1], followed by the remaining        */
/*    segments.  Two ints per segment.                                       */
/*  `segmentmarkerlist':  An array of segment markers; one int per segment.  */
/*                                                                           */
/*  `facetsegments':  The number of segments bounding each facet.  One int   */
/*    per facet.                                                             */
/*  `facetlist':  An array specifying which segments bound each facet.  The  */
/*    first facet's segments start at index [0].  The number of indices      */
/*    apportioned to each facet is the number of segments specified in the   */
/*    `facetsegments' array.                                                 */
/*  `facetmarkerlist':  An array of facet markers; one int per facet.        */
/*                                                                           */
/*  `holelist':  An array of holes.  The first hole's x, y, and z            */
/*    coordinates are at indices [0], [1], and [2], followed by the          */
/*    remaining holes.  Three REALs per hole.  Input only, although the      */
/*    pointer is copied to the output structure for your convenience.        */
/*                                                                           */
/*  `regionlist':  An array of regional attributes and volume constraints.   */
/*    The first constraint's x, y, and z coordinates are at indices [0],     */
/*    [1], and [2], followed by the regional attribute at index [3],         */
/*    followed by the maximum volume at index [4], followed by the remaining */
/*    regional attributes and volume constraints.  Five REALs per regional   */
/*    attribute/volume constraint.  Note that each regional attribute is     */
/*    used only if you select the `A' switch, and each volume constraint is  */
/*    used only if you select the `a' switch (with no number following), but */
/*    omitting one of these switches does not change the memory layout.      */
/*    Input only, although the pointer is copied to the output structure for */
/*    your convenience.                                                      */
/*                                                                           */
/*  `edgelist':  An array of edge endpoints.  The first edge's endpoints are */
/*    at indices [0] and [1], followed by the remaining edges.  Two ints per */
/*    edge.  Output only.                                                    */
/*  `edgemarkerlist':  An array of edge markers; one int per edge.  Output   */
/*    only.                                                                  */
/*  `normlist':  An array of normal vectors, used for infinite rays in       */
/*    Voronoi diagrams.  The first normal vector's x, y, and z magnitudes    */
/*    are at indices [0], [1], and [2], followed by the remaining vectors.   */
/*    For each finite edge in a Voronoi diagram, the normal vector written   */
/*    is the zero vector.  Three REALs per edge.  Output only.               */
/*                                                                           */
/*  `facelist':  An array of triangular face vertices.  The first face's     */
/*    vertices are at indices [0], [1], and [2], followed by the remaining   */
/*    faces.  Three ints per face.  Output only.                             */
/*  `facemarkerlist':  An array of face markers; one int per face.  Output   */
/*    only.                                                                  */
/*                                                                           */
/*                                                                           */
/*  Any input fields that Pyramid will examine must be initialized.          */
/*  Furthermore, for each output array that Pyramid will write to, you       */
/*  must either provide space by setting the appropriate pointer to point    */
/*  to the space you want the data written to, or you must initialize the    */
/*  pointer to NULL, which tells Pyramid to allocate space for the results.  */
/*  The latter option is preferable, because Pyramid always knows exactly    */
/*  how much space to allocate.  The former option is provided mainly for    */
/*  people who need to call Pyramid from Fortran code, though it also makes  */
/*  possible some nasty space-saving tricks, like writing the output to the  */
/*  same arrays as the input.                                                */
/*                                                                           */
/*  Pyramid will not free() any input or output arrays, including those it   */
/*  allocates itself; that's up to you.                                      */
/*                                                                           */
/*  Here's a guide to help you decide which fields you must initialize       */
/*  before you call tetrahedralize().                                        */
/*                                                                           */
/*  `in':                                                                    */
/*                                                                           */
/*    - `pointlist' must always point to a list of points; `numberofpoints'  */
/*      and `numberofpointattributes' must be properly set.                  */
/*      `pointmarkerlist' must either be set to NULL (in which case all      */
/*      markers default to zero), or must point to a list of markers.  If    */
/*      `numberofpointattributes' is not zero, `pointattributelist' must     */
/*      point to a list of point attributes.                                 */
/*    - If the `r' switch is used, `tetrahedronlist' must point to a list of */
/*      tetrahedra, and `numberoftetrahedra', `numberofcorners', and         */
/*      `numberoftetrahedronattributes' must be properly set.  If            */
/*      `numberoftetrahedronattributes' is not zero,                         */
/*      `tetrahedronattributelist' must point to a list of tetrahedron       */
/*      attributes.  If the `a' switch is used (with no number following),   */
/*      `tetrahedronvolumelist' must point to a list of tetrahedron volume   */
/*      constraints.  `neighborlist' may be ignored.                         */
/*    - If the `p' switch is used, `segmentlist' must point to a list of     */
/*      segment endpoints, `numberofsegments' must be properly set, and      */
/*      `segmentmarkerlist' must either be set to NULL (in which case all    */
/*      markers default to zero), or must point to a list of markers.        */
/*      `facetsegments' must point to a list of facet sizes, `facetlist'     */
/*      must point to a list of facet segments, `numberoffacets' must be     */
/*      properly set, and `facetmarkerlist' must either be set to NULL (in   */
/*      which case all markers default to zero), or must point to a list of  */
/*      markers.                                                             */
/*    - If the `p' switch is used without the `r' switch, then               */
/*      `numberofholes' and `numberofregions' must be properly set.  If      */
/*      `numberofholes' is not zero, `holelist' must point to a list of      */
/*      holes.  If `numberofregions' is not zero, `regionlist' must point to */
/*      a list of regional attributes/volume constraints.                    */
/*    - If the `p' switch is used, `holelist', `numberofholes',              */
/*      `regionlist', and `numberofregions' is copied to `out'.  (You can    */
/*      nonetheless get away with not initializing them if the `r' switch is */
/*      used.)                                                               */
/*    - `edgelist', `edgemarkerlist', `normlist', `numberofedges',           */
/*      `facelist', `facemarkerlist', and `numberoffaces' may be ignored.    */
/*                                                                           */
/*  `out':                                                                   */
/*                                                                           */
/*    - `pointlist' must be initialized (NULL or pointing to memory) unless  */
/*      the `N' switch is used.  `pointmarkerlist' must be initialized       */
/*      unless the `N' or `B' switch is used.  If `N' is not used and        */
/*      `in->numberofpointattributes' is not zero, `pointattributelist' must */
/*      be initialized.                                                      */
/*    - `tetrahedronlist' must be initialized unless the `E' switch is used. */
/*      `neighborlist' must be initialized if the `n' switch is used.  If    */
/*      the `E' switch is not used and (`in->numberofelementattributes' is   */
/*      not zero or the `A' switch is used), `elementattributelist' must be  */
/*      initialized.  `tetrahedronvolumelist' may be ignored.                */
/*    - `segmentlist' and `facetlist' must be initialized if the `p' or `c'  */
/*      switch is used, and the `P' switch is not used.  `segmentmarkerlist' */
/*      and `facetmarkerlist' must also be initialized under these           */
/*      circumstances unless the `B' switch is used.                         */
/*    - `edgelist' must be initialized if the `e' switch is used.            */
/*      `edgemarkerlist' must be initialized if the `e' switch is used and   */
/*      the `B' switch is not.                                               */
/*    - `facelist' must be initialized if the `f' switch is used.            */
/*      `facemarkerlist' must be initialized if the `f' switch is used and   */
/*      the `B' switch is not.                                               */
/*    - `holelist', `regionlist', `normlist', and all scalars may be ignored.*/
/*                                                                           */
/*  `vorout' (only needed if `v' switch is used):                            */
/*                                                                           */
/*    - `pointlist' must be initialized.  If `in->numberofpointattributes'   */
/*      is not zero, `pointattributelist' must be initialized.               */
/*      `pointmarkerlist' may be ignored.                                    */
/*    - `edgelist', `normlist', and `facelist' must all be initialized.      */
/*      `edgemarkerlist' and `facemarkerlist' may be ignored.                */
/*    - Everything else may be ignored.                                      */
/*                                                                           */
/*  After a call to tetrahedralize(), the valid fields of `out' and `vorout' */
/*  will depend, in an obvious way, on the choice of switches used.  Note    */
/*  that when the `p' switch is used, the pointers `holelist' and            */
/*  `regionlist' are copied from `in' to `out', but no new space is          */
/*  allocated; be careful that you don't free() the same array twice.  On    */
/*  the other hand, Pyramid will never copy the `pointlist' pointer (or any  */
/*  others); new space is allocated for `out->pointlist', or if the `N'      */
/*  switch is used, `out->pointlist' remains uninitialized.                  */
/*                                                                           */
/*  All of the meaningful `numberof' fields will be properly set; for        */
/*  instance, `numberofedges' will represent the number of edges in the      */
/*  triangulation whether or not the edges were written.  If a PLC is not    */
/*  used, `numberofsegments' will indicate the number of boundary edges, and */
/*  `numberoffacets' will indicate the number of boundary faces.             */
/*                                                                           */
/*****************************************************************************/

struct tetrahedralizeio {
  REAL *pointlist;                                               /* In / out */
  REAL *pointattributelist;                                      /* In / out */
  int *pointmarkerlist;                                          /* In / out */
  long numberofpoints;                                           /* In / out */
  int numberofpointattributes;                                   /* In / out */

  long *tetrahedronlist;                                         /* In / out */
  REAL *tetrahedronattributelist;                                /* In / out */
  REAL *tetrahedronvolumelist;                                    /* In only */
  long *neighborlist;                                            /* Out only */
  long numberoftetrahedra;                                       /* In / out */
  int numberofcorners;                                           /* In / out */
  int numberoftetrahedronattributes;                             /* In / out */

  long *segmentlist;                                             /* In / out */
  int *segmentmarkerlist;                                        /* In / out */
  long numberofsegments;                                         /* In / out */

  int *facetsegments;                                            /* In / out */
  long *facetlist;                                               /* In / out */
  int *facetmarkerlist;                                          /* In / out */
  long numberoffacets;                                           /* In / out */

  REAL *holelist;                        /* In / pointer to array copied out */
  long numberofholes;                                     /* In / copied out */

  REAL *regionlist;                      /* In / pointer to array copied out */
  long numberofregions;                                   /* In / copied out */

  long *edgelist;                                                /* Out only */
  int *edgemarkerlist;            /* Not used with Voronoi diagram; out only */
  REAL *normlist;                /* Used only with Voronoi diagram; out only */
  long numberofedges;                                            /* Out only */

  long *facelist;                                                /* Out only */
  int *facemarkerlist;            /* Not used with Voronoi diagram; out only */
  long numberoffaces;                                            /* Out only */
};

#ifdef ANSI_DECLARATORS
void tetrahedralize(char *, struct tetrahedralizeio *,
                    struct tetrahedralizeio *, struct tetrahedralizeio *);
#else /* not ANSI_DECLARATORS */
void tetrahedralize();
#endif /* not ANSI_DECLARATORS */
