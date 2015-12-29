/**
 * @author MrMaksymov / http://github.com/MrMaksymov
**/
//Can be used for straight extruding the shapes.
//shape2d is array of Vector2
//path2d is array of Vector2
//close bool= should this path be closed
THREE.StraightExtrudeGeometry = function ( shape2d, path2d, close )
{
    var scope = this;
    THREE.Geometry.call( this );
	this.type = 'StraightExtrudeGeometry';
	this.path2d = path2d;
	this.shape2d = shape2d;
    this.segments = [];
    var closed = this.path2d[0].clone().sub(this.path2d[path2d.length-1]).length()< 0.001;
    if(close && !closed){
        this.path2d.push(this.path2d[0]);
        closed = true;
    }

    if (this.path2d.length < 2) {console.warn("Can not extrude path < 2 points"); return;}
    if (this.shape2d.length < 3) {console.warn("Can not extrude shape < 3 points"); return;}


	for (var i=0; i <this.path2d.length-1; i++)
		this.segments.push(new THREE.Vector4(this.path2d[i].x,this.path2d[i].y,this.path2d[i+1].x,this.path2d[i+1].y));

	function LineXLine(data)
	{
		var line1 = data.line1;
		var line2 = data.line2;

		var A1 = line1.y - line1.w;
		var B1 = line1.z - line1.x;
		var C1 = line1.x * line1.w - line1.z * line1.y;
		var tSQR1 = Math.sqrt(A1 * A1 + B1 * B1);

		var A2 = line2.y - line2.w;
		var B2 = line2.z - line2.x;
		var C2 = line2.x * line2.w - line2.z * line2.y;
		var tSQR2 = Math.sqrt(A2 * A2 + B2 * B2);

        if (-A1/B1 == -A2/B2)
        {
            var lz = line2.x - line2.w - line2.y;
            var lw = line2.y + line2.z - line2.x;

            A2 = line2.y - lw;
            B2 = lz - line2.x;
            C2 = line2.x * lw - lz * line2.y;
            tSQR2 = Math.sqrt(A2 * A2 + B2 * B2);
        }

		var CP1 = C1 + data.delta1 * tSQR1;
		var CP2 = C2 + data.delta2 * tSQR2;

		return new THREE.Vector2((B1 * CP2 - B2 * CP1) / (A1 * B2 - A2 * B1) , (CP1 * A2 - CP2 * A1) / (A1 * B2 - A2 * B1));
	}

	for (var seg=0; seg <this.segments.length; seg++)
	{
		var c_seg = this.segments[seg];
		var c_ang = Math.atan2(c_seg.w - c_seg.y, c_seg.z - c_seg.x)+ Math.PI/2;

		var p_index = seg-1;
		if ((p_index == -1) && (closed)){ p_index = this.segments.length-1;}

		var n_index = seg+1;
		if ((n_index == this.segments.length) && (closed)){ n_index = 0;}

		for (var ps = 0; ps< this.shape2d.length-1; ps++)
		{
			var vp0 = new THREE.Vector3();
			var vp1 = new THREE.Vector3();
			{
				var locdata0= {"line1":c_seg , "line2": this.segments[p_index] , "delta1":this.shape2d[ps].x, "delta2":this.shape2d[ps].x};
				var locdata1= {"line1":c_seg , "line2": this.segments[p_index] , "delta1":this.shape2d[ps+1].x, "delta2":this.shape2d[ps+1].x};
				if (p_index == -1)
				{
					var pvec = new THREE.Vector4(c_seg.x,c_seg.y, c_seg.x + c_seg.y-c_seg.w, c_seg.y + c_seg.z - c_seg.x);
					locdata0.line2 = locdata1.line2 = pvec;
					locdata0.delta2 = locdata1.delta2 = 0;
				}
				var cross_locdata0 = LineXLine(locdata0);
                var cross_locdata1 = LineXLine(locdata1);
				vp0.set(cross_locdata0.x, this.shape2d[ps].y, cross_locdata0.y);
				vp1.set(cross_locdata1.x, this.shape2d[ps+1].y, cross_locdata1.y);
			}

			var vn0 = new THREE.Vector3();
			var vn1 = new THREE.Vector3();
			{
				var locdata0= {"line1":c_seg , "line2": this.segments[n_index] , "delta1":this.shape2d[ps].x, "delta2":this.shape2d[ps].x};
				var locdata1= {"line1":c_seg , "line2": this.segments[n_index] , "delta1":this.shape2d[ps+1].x, "delta2":this.shape2d[ps+1].x};
				if (n_index == this.segments.length)
				{
					var pvec = new THREE.Vector4(c_seg.z,c_seg.w, c_seg.z + c_seg.y - c_seg.w, c_seg.w + c_seg.z - c_seg.x);
					locdata0.line2 = locdata1.line2 = pvec;
					locdata0.delta2 = locdata1.delta2 = 0;
				}
                var cross_locdata0 = LineXLine(locdata0);
                var cross_locdata1 = LineXLine(locdata1);
				vn0.set(cross_locdata0.x, this.shape2d[ps].y, cross_locdata0.y);
				vn1.set(cross_locdata1.x, this.shape2d[ps+1].y, cross_locdata1.y);
			}

			var offset = this.vertices.length;
			this.vertices.push(vp0, vp1, vn0, vn1);
			this.faces.push(new THREE.Face3( offset+ 0, offset + 3, offset + 1 ) , new THREE.Face3( offset+ 0, offset + 2, offset + 3 ));
		}

	}

	if (!closed)
	{
		{
			var s_seg = this.segments[0];
			var s_ang = -Math.atan2(s_seg.w - s_seg.y, s_seg.z - s_seg.x)-Math.PI/2;
			var shape2dr = [];
			for (var i = 0; i< this.shape2d.length-1;i++)
				shape2dr.push(new THREE.Vector2(-this.shape2d[i].x, this.shape2d[i].y));
			var californiaShaped = new THREE.Shape( shape2dr );

			var geometrysg = new THREE.ShapeGeometry( californiaShaped );
			var mat = new THREE.Matrix4();
			var rotation = new THREE.Matrix4().makeRotationY(s_ang);
			var translation = new THREE.Matrix4().makeTranslation(s_seg.x, 0, s_seg.y);
			geometrysg.applyMatrix(rotation);
			geometrysg.applyMatrix(translation);
			this.merge(geometrysg, mat);
		}
		{
			var e_seg = this.segments[this.segments.length-1];
			var e_ang = -Math.atan2(e_seg.w - e_seg.y, e_seg.z - e_seg.x)+Math.PI/2;

			var californiaShaped = new THREE.Shape( this.shape2d );
			var geometryeg = new THREE.ShapeGeometry( californiaShaped );
			//
			var rotation2 = new THREE.Matrix4().makeRotationY(e_ang);
			var translation2 = new THREE.Matrix4().makeTranslation(e_seg.z, 0, e_seg.w);

			geometryeg.applyMatrix(rotation2);
			geometryeg.applyMatrix(translation2);
			this.merge(geometryeg, mat);
		}
	}

	//repaint UV
	this.faceVertexUvs[0] = [];
	var faces = this.faces;
	for (i = 0; i < this.faces.length ; i++)
	{
		var v1 = this.vertices[faces[i].a];
		var v2 = this.vertices[faces[i].b];
		var v3 = this.vertices[faces[i].c];

		var norm = new THREE.Triangle(v1,v2,v3).normal();
		norm.x = Math.abs(norm.x);norm.y = Math.abs(norm.y);norm.z = Math.abs(norm.z);
		var xplane = false;
		var yplane = false;
		var zplane = false;
		if ((norm.x <= norm.y) && (norm.z <= norm.y)) { xplane = true; yplane = false; zplane = true; }
		if ((norm.x <= norm.z) && (norm.y <= norm.z)) { xplane = true; yplane = true; zplane = false; }
		if ((norm.y <= norm.x) && (norm.z <= norm.x)) { xplane = false; yplane = true; zplane = true; }

		var mult=1;//если нужны мм в UV, то установить mult=1; 0.001 - это для теста, т.к. текстура 1х1 мм сливается)
		if ((xplane) && (zplane))
		{
			this.faceVertexUvs[0].push([new THREE.Vector2(v1.z*mult, v1.x*mult),new THREE.Vector2(v2.z*mult,v2.x*mult),new THREE.Vector2(v3.z*mult,v3.x*mult)]);
		}
		if ((yplane) && (zplane))
		{ //z y переставлены местами, т.к. текстура имеет вертикальное представление, а на плинтусе вероятно надо ее горизонтально рисовать
			this.faceVertexUvs[0].push([new THREE.Vector2(v1.y*mult, v1.z*mult), new THREE.Vector2(v2.y*mult,v2.z*mult),new THREE.Vector2(v3.y*mult,v3.z*mult)]);
		}
		if ((xplane) && (yplane))
		{
			this.faceVertexUvs[0].push([new THREE.Vector2(v1.y*mult, v1.x*mult), new THREE.Vector2(v2.y*mult,v2.x*mult),new THREE.Vector2(v3.y*mult,v3.x*mult)]);
		}
	}
	this.uvsNeedUpdate = true;
	this.computeFaceNormals();
	this.computeVertexNormals();    // requires correct face normals
    this.uvsNeedUpdate = true;
    this.verticesNeedUpdate = true;
    this.elementsNeedUpdate = true;
    this.computeCentroids;
    this.computeBoundingBox();
}


THREE.StraightExtrudeGeometry.prototype = Object.create( THREE.Geometry.prototype );
THREE.StraightExtrudeGeometry.prototype.constructor = THREE.StraightExtrudeGeometry;
THREE.StraightExtrudeGeometry.prototype.clone = function () {
	var parameters = this.parameters;
	return new THREE.ContentGeometry(
		parameters.shape2d,
		parameters.path2d
	);

};
