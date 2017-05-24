#pragma once

#include "sph_math.h"	

namespace SPH
{
	/** Fluid Point
	*/
	struct Point
	{
	public:
		fVector3		pos;					//λ��
		float			density;				//�ܶ�
		float			pressure;				//ѹ��
		fVector3		accel;					//���ٶ�
		float			temp;					//�¶�
		float			tempmar;				//�¶ȶ�ʱ��ĵ���
		float			thermal_conductivity;	//�ȴ���ϵ��

		fVector3		velocity;			
		fVector3		velocity_eval;		

		int				next;
	};

	/** point buffer
	*/
	class PointBuffer
	{
	public:
		void reset(unsigned int capcity);
		unsigned int size(void) const { return mFluidCounts; }
		unsigned int sizefluid(void) const { return mFluidPoints; }
		Point* get(unsigned int index) { return mFluidBuf+index; }
		const Point* get(unsigned int index) const { return mFluidBuf+index; }
		Point* AddPointReuse(void);
		Point* AddPointGhost(void);

	private:
		enum { ELEM_MAX=40960*2, };
		Point* mFluidBuf;
		unsigned int mFluidCounts;
		unsigned int mBufCapcity;
		unsigned int mFluidPoints;

	public:
		PointBuffer();
		virtual ~PointBuffer();
	};
}
