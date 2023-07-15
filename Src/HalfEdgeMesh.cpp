#include "HalfEdgeMesh.h"

using namespace Manifest_Math;

MFbool Manifest_Math::TrackEdge(HullHalfEdge* edge)
{
#ifdef DEBUG
	++edge->tail->referenceCount;
	DLOG(34, "adding edge: " << edge << " with vertex: " << edge->tail << "(" << edge->tail->vertex << ")  nReferences: " << edge->tail->referenceCount); 
#endif // DEBUG
	return true;
}
MFbool Manifest_Math::UntrackEdge(HullHalfEdge* edge)
{
#ifdef DEBUG
	--edge->tail->referenceCount;
	DLOG(34, "removing edge: " << edge << "(face: " << edge->face << ")  with vertex : " << edge->tail << " nReferences : " << edge->tail->referenceCount);
#endif // DEBUG
	return true;
}

void Manifest_Math::PrintFace(const HullFace* const face)
{
#ifdef DEBUG
	DLOG(43, "printing vertices of face: " << face);
	auto edge{ face->edge };
	do
	{
		DLOG(37, edge->tail->vertex);
		edge = edge->next;
	} while (edge != face->edge);
#endif // DEBUG
}

std::vector<HullHalfEdge*> Manifest_Math::GetFaceEdges(const HullFace* const face)
{
	DLOG(37, "Gathering edges of face: " << face);
	std::vector<HullHalfEdge*> result{ face->edge };
	for (auto edge{ face->edge->next }; edge != face->edge; edge = edge->next)
		result.emplace_back(edge);

	return result;
}

MFsize Manifest_Math::FaceVertexCount(const HullFace* const face)
{
	MFsize result{ 1 };
	for (auto edge{ face->edge->next }; edge != face->edge; edge = edge->next)
		++result;

	return result;
}

MFplane Manifest_Math::CalculateNormalizedFacePlane(const HullFace* const face)
{
	MFplane result;
	const auto planarVertexCount{ FaceVertexCount(face) };
	if (planarVertexCount > 3)
		result = NewellPlane(planarVertexCount, face);
	else
	{
		assert(planarVertexCount == 3);
		MFtriangle triangle;
		triangle.vertices[0] = face->edge->tail->vertex;
		triangle.vertices[1] = face->edge->next->tail->vertex;
		triangle.vertices[2] = face->edge->prev->tail->vertex;
		result = Normalize(CalculateSurfacePlane(triangle));
	}

	return result;
}

MFplane Manifest_Math::NewellPlane(const MFu32& planarEdgeCount, const HullFace* const planarFace)
{
	MFvec3 normal{ 0 };
	auto planarEdge{ planarFace->edge };
	for (auto edge{ 0 }; edge < planarEdgeCount; ++edge)
	{
		const auto& currentVertex{ planarEdge->tail->vertex };
		const auto& nextVertex{ planarEdge->next->tail->vertex };
		normal.x += (currentVertex.y - nextVertex.y) * (currentVertex.z + nextVertex.z);
		normal.y += (currentVertex.z - nextVertex.z) * (currentVertex.x + nextVertex.x);
		normal.z += (currentVertex.x - nextVertex.x) * (currentVertex.y + nextVertex.y);
		planarEdge = planarEdge->next;
	}
	normal = Normalize(normal);
	const auto offset{ Dot(-normal,planarFace->edge->tail->vertex) };
	return MFplane{ normal,offset };
}

//uses Eurler's Formula to allocate double the amount required baseline
void HalfEdgeMesh::AllocateBuffers(const MFsize& numberVertices)
{
	vertexBuffer = new HullVertex[numberVertices * 2];
	edgeBuffer = new HullHalfEdge[(3 * numberVertices - 6) * 2 * 2];//half edges,2*2
	faceBuffer = new HullFace[(2 * numberVertices - 4) * 2];
}

void HalfEdgeMesh::DeallocateBuffers()
{
	if (vertexBuffer)
	{
		delete[] vertexBuffer;
		vertexBuffer = nullptr;
	}
	if (edgeBuffer)
	{
		delete[] edgeBuffer;
		edgeBuffer = nullptr;
	}
	if (faceBuffer)
	{
		delete[] faceBuffer;
		faceBuffer = nullptr;
	}
}

HullVertex* HalfEdgeMesh::AllocateVertex()
{
	HullVertex* result;
	if (vertexFreeList.size())
	{
		result = vertexFreeList.top();
		vertexFreeList.pop();
	}
	else
		result = &vertexBuffer[allocatedVertices++];

	return result;
}

HullHalfEdge* HalfEdgeMesh::AllocateHalfEdge()
{
	HullHalfEdge* result{ nullptr };
	//change this in a second to test if you can do a boolean set test
	if (edgeFreeList.size())
	{
		result = edgeFreeList.top();
		edgeFreeList.pop();
		DLOG(35, "Allocating reuse edge: " << result);
	}
	else
	{
		result = &edgeBuffer[allocatedEdges++];
		DLOG(36, "Allocating new edge: " << result);
	}
	++edgeCount;
	return new(result)HullHalfEdge;
}

HullFace* HalfEdgeMesh::AllocateFace()
{
	HullFace* result{ nullptr };
	if (faceFreeList.size())
	{
		result = faceFreeList.top();
		faceFreeList.pop();
		DLOG(37, "Allocating reuse face: " << result);
	}
	else
	{
		result = &faceBuffer[allocatedFaces++];
		DLOG(33, "Allocating new face: " << result);
	}

	return new(result)HullFace;
}
void HalfEdgeMesh::DeallocateVertex(HullVertex* vertex)
{
	DLOG(44, "removing vertex: " << vertex << " " << vertex->vertex);
	//vertex should only be removed when no edges use it
	assert(vertex->referenceCount == 0);
	//only needs to run if vertex was originally added to hull
	if (vertex->prev)//will always have prev/next pointers if so
		RemoveHalfEdgeMeshFeature<HullVertex>(vertexCount, vertices, vertex);
	vertexFreeList.push(vertex);
}
void HalfEdgeMesh::DeallocateHalfEdge(HullHalfEdge* edge)
{
	UntrackEdge(edge);
	edge->deallocated = true;
	edgeFreeList.push(edge);
	--edgeCount;
}
void HalfEdgeMesh::DeallocateFace(HullFace* face)
{
	DLOG(34, "Deallocating face: " << face);
	RemoveHalfEdgeMeshFeature<HullFace>(faceCount, faces, face);
	face->deallocated = true;
	faceFreeList.push(face);
}
