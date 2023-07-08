#pragma once
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <list>
#include <type_traits>

#include "Plane.h"
#include "Geometry.h"

namespace Manifest_Math
{
	struct HullVertex
	{
		HullVertex* prev{ nullptr };
		HullVertex* next{ nullptr };

		MFpoint3 vertex{ 0 };
		MFu32 referenceCount{ 0 };
	};
	using ConflictList = std::vector<HullVertex>;

	struct HullFace;

	struct HullHalfEdge
	{
		HullHalfEdge* prev{ nullptr };
		HullHalfEdge* next{ nullptr };
		HullHalfEdge* twin{ nullptr };
		HullVertex* tail{ nullptr };
		HullFace* face;

		MFbool deallocated{ false };
	};
	//TTY tracking functions - debug mode	
	MFbool TrackEdge(HullHalfEdge* edge);
	MFbool UntrackEdge(HullHalfEdge* edge);

	struct HullFace
	{
		HullFace* prev{ nullptr };
		HullFace* next{ nullptr };
		HullHalfEdge* edge{ nullptr };

		std::list<HullVertex*> conflistList;
		MFbool deallocated{ false };
	};
	void PrintFace(const HullFace* const face);
	std::vector<HullHalfEdge*> GetFaceEdges(const HullFace* const face);
	MFsize FaceVertexCount(const HullFace* const face);
	MFplane CalculateNormalizedFacePlane(const HullFace* const face);
	MFplane NewellPlane(const MFu32& planarEdgeCount, const HullFace* const planarFace);

	class HalfEdgeMesh
	{	
		public:
			//allocates inital data stores
			void AllocateBuffers(const MFsize& numberVertices);
			void DeallocateBuffers();
			//de/allocates hull structures as required		
			void DeallocateVertex(HullVertex* vertex);
			void DeallocateHalfEdge(HullHalfEdge* edge);
			void DeallocateFace(HullFace* face);
			HullVertex* AllocateVertex();
			HullHalfEdge* AllocateHalfEdge();
			HullFace* AllocateFace();

			HullVertex* vertices{ nullptr };			
			HullFace* faces{ nullptr };
			MFsize vertexCount;
			MFsize edgeCount;
			MFsize faceCount;
			MFfloat CONVEXITY_EPSILON{ FLT_EPSILON };				
		private:						
			//allocation buffers
			HullVertex* vertexBuffer{ nullptr };
			HullHalfEdge* edgeBuffer{ nullptr };
			HullFace* faceBuffer{ nullptr };
			//represents total number of allocations, only increases
			//once allocated either on hull or in free list
			MFsize allocatedVertices{ 0 };
			MFsize allocatedEdges{ 0 };
			MFsize allocatedFaces{ 0 };
			//allocates off freelist first if possible
			std::stack<HullVertex*> vertexFreeList;
			std::stack<HullHalfEdge*> edgeFreeList;
			std::stack<HullFace*> faceFreeList;			
	};
	template<typename T>
	void AddHalfEdgeMeshFeature(MFsize& featureCount, T*& meshFeatures, T* additionalFeature)
	{		
		if (meshFeatures)
		{
			auto back{ meshFeatures->prev };
			back->next = additionalFeature;
			additionalFeature->prev = back;			
		}
		else		
			meshFeatures = additionalFeature;	
		
		additionalFeature->next = meshFeatures;
		meshFeatures->prev = additionalFeature;
		
		++featureCount;
	}
	template<typename T>
	void RemoveHalfEdgeMeshFeature(MFsize& featureCount, T*& meshFeatures, T* removedFeature)
	{				
		if (meshFeatures == removedFeature)
			meshFeatures = meshFeatures->next;
		
		removedFeature->prev->next = removedFeature->next;
		removedFeature->next->prev = removedFeature->prev;			

		--featureCount;
	}
}