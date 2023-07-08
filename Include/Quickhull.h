#pragma once
#include "HalfEdgeMesh.h"
#include "Geometry.h"
#include "Simplex.h"
#include "Vector2.h"

namespace Manifest_Math
{
	const MFfloat EPSILON_3D(const std::vector<MFpoint3>& pointCloud);

	HalfEdgeMesh QuickHull(std::vector<MFpoint3> pointCloud);
	//performs complete convexity check
	MFbool AssertConvexity(const HalfEdgeMesh& mesh);	

	//construction functions
	Simplex_T<MFpoint3> ConstructInitialSimplex(const std::vector<MFpoint3>& vertices, MFvec4& CONTRUCTION_NORMAL_OUT);
	void ConstructInitialHull(const Simplex_T<MFpoint3>& simplex, const MFvec4& CONTRUCTION_NORMAL_IN, std::vector<MFpoint3>& pointCloud, HalfEdgeMesh& mesh);
	void AddPointToHull(HalfEdgeMesh& mesh, const MFpoint3& newPoint);
	HullVertex* NextConflictVertex(HalfEdgeMesh& mesh, HullFace*& conflictFace);
	void AddVertexToHull(HalfEdgeMesh& mesh, HullFace* conflictFace, HullVertex* vertex);
	std::vector<HullFace*> BuildHorizon(HalfEdgeMesh& mesh, HullVertex* const vertex, std::list<HullHalfEdge*>& horizon, HullFace* conflictFace);
	void HorizonDFS(HullFace* face, const HullVertex* const vertex, std::unordered_set<HullFace*>& visitedFaces, std::list<HullHalfEdge*>& horizon, std::vector<HullFace*>& visibleFaces, const MFfloat epsilon);
	void BuildNewFaces(HalfEdgeMesh& mesh, HullVertex* const vertex, std::vector<HullFace*>& newFaces, std::list<HullHalfEdge*>& horizon);
	void MergeFaces(HalfEdgeMesh& mesh, std::vector<HullFace*>& newFaces, HullFace* conflictFace);
	MFbool CheckFaceConvexity(const HullFace* const centroidFace, const HullFace* const planarFace, const MFfloat& HULL_EPSILON);
	void UpdateHullFaces(HalfEdgeMesh& mesh, std::vector<HullFace*>& newFaces, std::vector<HullFace*>& visibleFaces, HullFace* conflictFace);
	void ResolveOrphans(HalfEdgeMesh& mesh, std::vector<HullFace*>& newFaces, HullFace* conflictFace);
};