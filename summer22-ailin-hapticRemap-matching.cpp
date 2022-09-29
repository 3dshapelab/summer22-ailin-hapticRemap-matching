// this script aims to use probe adjustment task to measure the gain/strength of either texture or disparity as a depth cue
#include "stdafx.h"

void initPreviewStimulus(double textDepth, double dispDepth) {

	if(abs(textDepth - dispDepth) < 0.1){
		display_distance_jittered_std = display_distance + jitter_z_std;
		buildVertices_congruent(true, textDepth, normalizer_to_uv_std);
	}else{	
		display_distance_jittered_std = display_distance + jitter_z_std;
		buildVertices_incongruent(true, textDepth, dispDepth, display_distance_jittered_std, normalizer_to_uv_std);
	}

}


void buildVertices_congruent(bool isStandard, double shapeDepth, double textNormalizer) {

	std::vector<GLfloat> vertices_vec;
	std::vector<GLfloat> texcoors_vec;
	std::vector<GLfloat> colors_vec;
	std::vector<GLfloat> normals_vec;
	std::vector<GLuint> indices_draw_triangle_vec;

	std::vector<Vector3d> vertContainer_Rcontour_Clpeye;
	std::vector<Vector3d> vertContainer_Lcontour_Clpeye;

	
	double step_size_width = stimulus_width / (double)(nr_points - 1);
	double step_size_height = stimulus_height / (double)(nr_points - 1);

	GLuint i_ind = 0;

	double x, y, z, y_prev, z_prev, u, v;
	y_prev = -stimulus_height / 2; 
	z_prev = 0; 
	double total_distance_y = 0; //tracks the distance along y/z axis, approximate the "diameter" of the ellipse

	double normal_x, normal_y, normal_z;

	double x_Rcontour = stimulus_visiblewidth / 2;
	double x_Lcontour = -stimulus_visiblewidth / 2;
	

	for (int j = 0; j < nr_points; j++) {  // 

		y = -stimulus_height / 2 + j * step_size_height;
		z = shapeDepth * cos(M_PI * y / stimulus_height);

		total_distance_y = total_distance_y + sqrt(pow(y - y_prev, 2) + pow(z - z_prev, 2));
		v = total_distance_y / textNormalizer + v_offset; //v coordinate

		// shading should be consistent with texture depth
		normal_x = 0;
		normal_y = shapeDepth * sin(M_PI * y / stimulus_height) * M_PI / stimulus_height;
		normal_z = 1;

		vertContainer_Rcontour_Clpeye.push_back(Vector3d(x_Rcontour, y, z)); 
		vertContainer_Lcontour_Clpeye.push_back(Vector3d(x_Lcontour, y, z)); 

		for (int i = 0; i < nr_points; i++) { //

			x = -stimulus_width / 2 + i * step_size_width;
			u = (x + stimulus_width / 2) / textNormalizer + u_offset; //u coordinate. 

			//step 1: build the meshgrid using vertices, 
			vertices_vec.push_back(x);
			vertices_vec.push_back(y);
			vertices_vec.push_back(z);

			colors_vec.push_back(1);
			colors_vec.push_back(0);
			colors_vec.push_back(0);

			texcoors_vec.push_back(u);
			texcoors_vec.push_back(v);

			normals_vec.push_back(normal_x );
			normals_vec.push_back(normal_y );
			normals_vec.push_back(normal_z );

			//step 2: create an array/vector that store how the triangles should be drawn

			// construct the triangle indices to be drawn
			if (i < nr_points - 1 && j < nr_points - 1) {

				indices_draw_triangle_vec.push_back(i_ind);
				indices_draw_triangle_vec.push_back(i_ind + 1);
				indices_draw_triangle_vec.push_back(i_ind + nr_points);

				indices_draw_triangle_vec.push_back(i_ind + nr_points);
				indices_draw_triangle_vec.push_back(i_ind + 1);
				indices_draw_triangle_vec.push_back(i_ind + nr_points + 1);
				//ind = ind + 6;
			}

			i_ind++;
		}

		y_prev = y; z_prev = z;
	}

	if(isStandard){///////////////////////// building standard /////////////////////////	
		
		vertices_vec_std.clear();
		colors_vec_std.clear();
		texcoors_vec_std.clear();
		indices_draw_triangle_vec_std.clear();
		normals_vec_std.clear();

		vertices_vec_std = vertices_vec;
		colors_vec_std = colors_vec;
		texcoors_vec_std = texcoors_vec;
		indices_draw_triangle_vec_std = indices_draw_triangle_vec;
		normals_vec_std = normals_vec;

		vertContainer_std_Rcontour.clear();
		vertContainer_std_Lcontour.clear();

		vertContainer_std_Rcontour = vertContainer_Rcontour_Clpeye;
		vertContainer_std_Lcontour = vertContainer_Lcontour_Clpeye;

	}else{

		vertices_vec_cmp.clear();
		colors_vec_cmp.clear();
		texcoors_vec_cmp.clear();
		indices_draw_triangle_vec_cmp.clear();
		normals_vec_cmp.clear();

		vertices_vec_cmp = vertices_vec;
		colors_vec_cmp = colors_vec;
		texcoors_vec_cmp = texcoors_vec;
		indices_draw_triangle_vec_cmp = indices_draw_triangle_vec;
		normals_vec_cmp = normals_vec;


		vertContainer_cmp_Rcontour.clear();
		vertContainer_cmp_Lcontour.clear();

		vertContainer_cmp_Rcontour = vertContainer_Rcontour_Clpeye;
		vertContainer_cmp_Lcontour = vertContainer_Lcontour_Clpeye;
	
	}

}


void buildVertices_incongruent(bool isStandard, double textDepth, double dispDepth, double displayDist, double textNormalizer) {
	
	double step_size_width = stimulus_width / (double)(nr_points - 1);
	double step_size_height = stimulus_height / (double)(nr_points - 1);

	GLuint i_ind = 0;

	std::vector<GLfloat> vertices_vec;
	std::vector<GLfloat> texcoors_vec;
	std::vector<GLfloat> colors_vec;
	std::vector<GLfloat> normals_vec;
	std::vector<GLuint> indices_draw_triangle_vec;

	std::vector<Vector3d> vertContainer_Rcontour_Clpeye;
	std::vector<Vector3d> vertContainer_Lcontour_Clpeye;

	double x_t, y_t, z_t, y_t_prev, z_t_prev, u, v;
	double w, x_d, y_d, z_d;
	y_t_prev = -stimulus_height / 2; 
	z_t_prev = 0; 
	double total_distance_y_t = 0; //tracks the distance along y_t/z_t axis, approximate the "diameter" of the ellipse
	
	double c_cos = M_PI;
	double l = -(displayDist - dispDepth);

	double x_t_Rcontour = stimulus_visiblewidth / 2;
	double x_t_Lcontour = -stimulus_visiblewidth / 2;
	double x_d_Rcontour_Leye, x_d_Rcontour_Reye, x_d_Lcontour_Leye, x_d_Lcontour_Reye;
	double x_d_Rcontour_Clpeye, x_d_Lcontour_Clpeye;

	double normal_x, normal_y, normal_z, normal_norm;
	
	for (int j = 0; j < nr_points; j++) {  // 

		y_t = -stimulus_height / 2 + j * step_size_height;
		z_t = textDepth * cos(M_PI * y_t / stimulus_height);

		// shading should be consistent with texture depth
		normal_x = 0;
		normal_y = textDepth * sin(M_PI * y_t / stimulus_height) * M_PI / stimulus_height;
		normal_z = 1;
		//normal_norm = sqrt(pow(normal_y,2) + pow(normal_z, 2));

		total_distance_y_t = total_distance_y_t + sqrt(pow(y_t - y_t_prev, 2) + pow(z_t - z_t_prev, 2));
		v = total_distance_y_t / textNormalizer + v_offset; //v coordinate
				
		if(abs(y_t) < 0.01){
			y_d = y_t;
			z_d = dispDepth * cos(M_PI * y_d / stimulus_height) ;

		}else if(abs( abs(y_t)- stimulus_height / 2) < 0.01){
			y_d = y_t;
			z_d = z_t;

		}else{
			c_cos = M_PI * y_t /(stimulus_height * (z_t - l));
			z_d = NewtonSolver_Cosine(stimulus_height, dispDepth, c_cos, l, y_t, z_t);
			w = (z_d - l)/(z_t - l);
			y_d = w * y_t;

		}

		w = (z_d - l)/(z_t - l);

		//x_d_Rcontour_Leye = -(interoculardistance / 2) + w * (x_t_Rcontour + interoculardistance / 2);
		//x_d_Lcontour_Leye = -(interoculardistance / 2) + w * (x_t_Lcontour + interoculardistance / 2);
		//x_d_Rcontour_Reye = (interoculardistance / 2) + w * (x_t_Rcontour - interoculardistance / 2);
		//x_d_Lcontour_Reye = (interoculardistance / 2) + w * (x_t_Lcontour - interoculardistance / 2);
		//vertContainer_std_Rcontour_Leye.push_back(Vector3d(x_d_Rcontour_Leye, y_d, z_d)); 
		//vertContainer_std_Rcontour_Reye.push_back(Vector3d(x_d_Rcontour_Reye, y_d, z_d));
		//vertContainer_std_Lcontour_Leye.push_back(Vector3d(x_d_Lcontour_Leye, y_d, z_d)); 
		//vertContainer_std_Lcontour_Reye.push_back(Vector3d(x_d_Lcontour_Reye, y_d, z_d));
		
		x_d_Rcontour_Clpeye = w * x_t_Rcontour;
		x_d_Lcontour_Clpeye = w * x_t_Lcontour ;

		vertContainer_Rcontour_Clpeye.push_back(Vector3d(x_d_Rcontour_Clpeye, y_d, z_d)); 
		vertContainer_Lcontour_Clpeye.push_back(Vector3d(x_d_Lcontour_Clpeye, y_d, z_d));

		for (int i = 0; i < nr_points; i++) { //

			x_t = -stimulus_width / 2 + i * step_size_width;
			u = (x_t + stimulus_width / 2) / textNormalizer + u_offset; //u coordinate. 

			x_d = x_t; 
	
			//step 1: build the meshgrid using vertices, 

			vertices_vec.push_back(x_d);
			vertices_vec.push_back(y_d);
			vertices_vec.push_back(z_d);

			colors_vec.push_back(1);
			colors_vec.push_back(0);
			colors_vec.push_back(0);

			texcoors_vec.push_back(u);
			texcoors_vec.push_back(v);

			normals_vec.push_back(normal_x);
			normals_vec.push_back(normal_y);
			normals_vec.push_back(normal_z);

			//step 2: create an array/vector that store how the triangles should be drawn

			// construct the triangle indices to be drawn
			if (i < nr_points - 1 && j < nr_points - 1) {

				indices_draw_triangle_vec.push_back(i_ind);
				indices_draw_triangle_vec.push_back(i_ind + 1);
				indices_draw_triangle_vec.push_back(i_ind + nr_points);

				indices_draw_triangle_vec.push_back(i_ind + nr_points);
				indices_draw_triangle_vec.push_back(i_ind + 1);
				indices_draw_triangle_vec.push_back(i_ind + nr_points + 1);
				//ind = ind + 6;
			}

			i_ind++;

		}
				
		y_t_prev = y_t; z_t_prev = z_t;
	}

	if(isStandard){///////////////////////// building standard /////////////////////////	
		
		vertices_vec_std.clear();
		colors_vec_std.clear();
		texcoors_vec_std.clear();
		indices_draw_triangle_vec_std.clear();
		normals_vec_std.clear();

		vertices_vec_std = vertices_vec;
		colors_vec_std = colors_vec;
		texcoors_vec_std = texcoors_vec;
		indices_draw_triangle_vec_std = indices_draw_triangle_vec;
		normals_vec_std = normals_vec;


		vertContainer_std_Rcontour.clear();
		vertContainer_std_Lcontour.clear();

		vertContainer_std_Rcontour = vertContainer_Rcontour_Clpeye;
		vertContainer_std_Lcontour = vertContainer_Lcontour_Clpeye;


	}else{

		vertices_vec_cmp.clear();
		colors_vec_cmp.clear();
		texcoors_vec_cmp.clear();
		indices_draw_triangle_vec_cmp.clear();
		normals_vec_cmp.clear();

		vertices_vec_cmp = vertices_vec;
		colors_vec_cmp = colors_vec;
		texcoors_vec_cmp = texcoors_vec;
		indices_draw_triangle_vec_cmp = indices_draw_triangle_vec;
		normals_vec_cmp = normals_vec;


		vertContainer_cmp_Rcontour.clear();
		vertContainer_cmp_Lcontour.clear();

		vertContainer_cmp_Rcontour = vertContainer_Rcontour_Clpeye;
		vertContainer_cmp_Lcontour = vertContainer_Lcontour_Clpeye;
	
	}

}


double NewtonSolver_fz(double var_z, double D, double C, double l){
	double val = var_z/D - cos(C*(var_z - l));
	return val;
}

double NewtonSolver_dfz(double var_z, double D, double C, double l){
	double val = 1/D + sin(C*(var_z - l)) * C;
	return val;
}

double NewtonSolver_Cosine(double theHeight, double newDeth_zDenom, double constCos, double l_translatedDist, double y0, double z0){

	double z_new, z, f_z, df_z;
	double C_cos = M_PI * y0/(theHeight * (z0 - l_translatedDist));

	z = z0;

	for (int i = 0; i < 100; i++){

		f_z = NewtonSolver_fz(z, newDeth_zDenom, C_cos, l_translatedDist);
		df_z = NewtonSolver_dfz(z, newDeth_zDenom, C_cos, l_translatedDist);

		if(abs(f_z) < 1e-10){
			//cout << "iteration: " << i << "      success!  fz = " << f_z << endl; 
			break;
		}else if(abs(df_z) < 1e-10){
			//cout << "iteration: " << i << "      uhhhhh  fz = " << f_z << endl; 
			break;
		}else{		
			z_new = z - f_z/df_z;
			z = z_new;
		}
	}

	if(abs(z-z0) > 40)
		z = 0;

	return z;

}

void drawPanels(bool isStandard, double displayDist, double dispDepth){

	int n;
	float panel_width = 40;
	float panel_height_extra = 20;

	glPushMatrix();
	glLoadIdentity();
	glTranslated(0, 0, displayDist - dispDepth + 2);

	if(isStandard){

		n = int(vertContainer_std_Rcontour.size());

		if(n > 0){

			// Right panels
			glBegin(GL_QUAD_STRIP);

			glVertex3f(vertContainer_std_Rcontour.at(0)[0] + panel_width,		vertContainer_std_Rcontour.at(0)[1] - panel_height_extra,			vertContainer_std_Rcontour.at(0)[2]); //0
			glVertex3f(vertContainer_std_Rcontour.at(0)[0],					vertContainer_std_Rcontour.at(0)[1] - panel_height_extra,			vertContainer_std_Rcontour.at(0)[2]); //1

			for (int i = 0; i < n; i++)
			{	
				glVertex3f(vertContainer_std_Rcontour.at(i)[0] + panel_width,		vertContainer_std_Rcontour.at(i)[1],			vertContainer_std_Rcontour.at(i)[2]); //0
				glVertex3f(vertContainer_std_Rcontour.at(i)[0],					vertContainer_std_Rcontour.at(i)[1],			vertContainer_std_Rcontour.at(i)[2]); //1

			}	

			glVertex3f(vertContainer_std_Rcontour.at(n-1)[0] + panel_width,		vertContainer_std_Rcontour.at(n-1)[1] + panel_height_extra,			vertContainer_std_Rcontour.at(n-1)[2]); //0
			glVertex3f(vertContainer_std_Rcontour.at(n-1)[0],					vertContainer_std_Rcontour.at(n-1)[1] + panel_height_extra,			vertContainer_std_Rcontour.at(n-1)[2]); //1

			glEnd();

			// Left panels
			glBegin(GL_QUAD_STRIP);

			glVertex3f(vertContainer_std_Lcontour.at(0)[0],		vertContainer_std_Lcontour.at(0)[1] - panel_height_extra,			vertContainer_std_Lcontour.at(0)[2]); //0
			glVertex3f(vertContainer_std_Lcontour.at(0)[0] - panel_width,					vertContainer_std_Lcontour.at(0)[1] - panel_height_extra,			vertContainer_std_Lcontour.at(0)[2]); //1

			for (int i = 0; i < n; i++)
			{	
				glVertex3f(vertContainer_std_Lcontour.at(i)[0],		vertContainer_std_Lcontour.at(i)[1],			vertContainer_std_Lcontour.at(i)[2]); //0
				glVertex3f(vertContainer_std_Lcontour.at(i)[0] - panel_width,					vertContainer_std_Lcontour.at(i)[1],			vertContainer_std_Lcontour.at(i)[2]); //1

			}	

			glVertex3f(vertContainer_std_Lcontour.at(n-1)[0],		vertContainer_std_Lcontour.at(n-1)[1] + panel_height_extra,			vertContainer_std_Lcontour.at(n-1)[2]); //0
			glVertex3f(vertContainer_std_Lcontour.at(n-1)[0] - panel_width,					vertContainer_std_Lcontour.at(n-1)[1] + panel_height_extra,			vertContainer_std_Lcontour.at(n-1)[2]); //1

			glEnd();
		}

	}else{

		n = int(vertContainer_cmp_Rcontour.size());

		if(n > 0){
			// Right panels
			glBegin(GL_QUAD_STRIP);

			glVertex3f(vertContainer_cmp_Rcontour.at(0)[0] + panel_width,		vertContainer_cmp_Rcontour.at(0)[1] - panel_height_extra,			vertContainer_cmp_Rcontour.at(0)[2]); //0
			glVertex3f(vertContainer_cmp_Rcontour.at(0)[0],					vertContainer_cmp_Rcontour.at(0)[1] - panel_height_extra,			vertContainer_cmp_Rcontour.at(0)[2]); //1
			for (int i = 0; i < n; i++)
			{	
				glVertex3f(vertContainer_cmp_Rcontour.at(i)[0] + panel_width,		vertContainer_cmp_Rcontour.at(i)[1],			vertContainer_cmp_Rcontour.at(i)[2]); 
				glVertex3f(vertContainer_cmp_Rcontour.at(i)[0],					vertContainer_cmp_Rcontour.at(i)[1],			vertContainer_cmp_Rcontour.at(i)[2]); 

			}	
			glVertex3f(vertContainer_cmp_Rcontour.at(n-1)[0] + panel_width,		vertContainer_cmp_Rcontour.at(n-1)[1] + panel_height_extra,			vertContainer_cmp_Rcontour.at(n-1)[2]); //0
			glVertex3f(vertContainer_cmp_Rcontour.at(n-1)[0],					vertContainer_cmp_Rcontour.at(n-1)[1] + panel_height_extra,			vertContainer_cmp_Rcontour.at(n-1)[2]); //1

			glEnd();

			// Left panels
			glBegin(GL_QUAD_STRIP);

			glVertex3f(vertContainer_cmp_Lcontour.at(0)[0],		vertContainer_cmp_Lcontour.at(0)[1] - panel_height_extra,			vertContainer_cmp_Lcontour.at(0)[2]); //0
			glVertex3f(vertContainer_cmp_Lcontour.at(0)[0] - panel_width,					vertContainer_cmp_Lcontour.at(0)[1] - panel_height_extra,			vertContainer_cmp_Lcontour.at(0)[2]); //1
			for (int i = 0; i < n; i++)
			{	
				glVertex3f(vertContainer_cmp_Lcontour.at(i)[0],		vertContainer_cmp_Lcontour.at(i)[1],			vertContainer_cmp_Lcontour.at(i)[2]); //0
				glVertex3f(vertContainer_cmp_Lcontour.at(i)[0] - panel_width,					vertContainer_cmp_Lcontour.at(i)[1],			vertContainer_cmp_Lcontour.at(i)[2]); //1

			}	
			glVertex3f(vertContainer_cmp_Lcontour.at(n-1)[0],		vertContainer_cmp_Lcontour.at(n-1)[1] + panel_height_extra,			vertContainer_cmp_Lcontour.at(n-1)[2]); //0
			glVertex3f(vertContainer_cmp_Lcontour.at(n-1)[0] - panel_width,					vertContainer_cmp_Lcontour.at(n-1)[1] + panel_height_extra,			vertContainer_cmp_Lcontour.at(n-1)[2]); //1

			glEnd();			

		}		
	
	}
	glPopMatrix();
}

void drawPanel_Leye(bool isStandard, double displayDist, double dispDepth){

	float panel_width = 40;
	float panel_height_extra = 20;
	int n = int(vertContainer_std_Rcontour_Leye.size());
	
	if(n > 0){

		glPushMatrix();	
		glLoadIdentity();
		glTranslated(0, 0, displayDist - dispDepth + 2);
		glColor3f(0.0f, 0.0f, 0.0f);

		glBegin(GL_QUAD_STRIP);

		glVertex3f(vertContainer_std_Rcontour_Leye.at(0)[0] + panel_width,		vertContainer_std_Rcontour_Leye.at(0)[1] - panel_height_extra,			vertContainer_std_Rcontour_Leye.at(0)[2]); //0
		glVertex3f(vertContainer_std_Rcontour_Leye.at(0)[0],					vertContainer_std_Rcontour_Leye.at(0)[1] - panel_height_extra,			vertContainer_std_Rcontour_Leye.at(0)[2]); //1

		for (int i = 0; i < n; i++)
		{	
			glVertex3f(vertContainer_std_Rcontour_Leye.at(i)[0] + panel_width,		vertContainer_std_Rcontour_Leye.at(i)[1],			vertContainer_std_Rcontour_Leye.at(i)[2]); //0
			glVertex3f(vertContainer_std_Rcontour_Leye.at(i)[0],					vertContainer_std_Rcontour_Leye.at(i)[1],			vertContainer_std_Rcontour_Leye.at(i)[2]); //1

		}	

		glVertex3f(vertContainer_std_Rcontour_Leye.at(n-1)[0] + panel_width,		vertContainer_std_Rcontour_Leye.at(n-1)[1] + panel_height_extra,			vertContainer_std_Rcontour_Leye.at(n-1)[2]); //0
		glVertex3f(vertContainer_std_Rcontour_Leye.at(n-1)[0],					vertContainer_std_Rcontour_Leye.at(n-1)[1] + panel_height_extra,			vertContainer_std_Rcontour_Leye.at(n-1)[2]); //1

		glEnd();


		glBegin(GL_QUAD_STRIP);

		glVertex3f(vertContainer_std_Lcontour_Leye.at(0)[0],		vertContainer_std_Lcontour_Leye.at(0)[1] - panel_height_extra,			vertContainer_std_Lcontour_Leye.at(0)[2]); //0
		glVertex3f(vertContainer_std_Lcontour_Leye.at(0)[0] - panel_width,					vertContainer_std_Lcontour_Leye.at(0)[1] - panel_height_extra,			vertContainer_std_Lcontour_Leye.at(0)[2]); //1

		for (int i = 0; i < n; i++)
		{	
			glVertex3f(vertContainer_std_Lcontour_Leye.at(i)[0],		vertContainer_std_Lcontour_Leye.at(i)[1],			vertContainer_std_Lcontour_Leye.at(i)[2]); //0
			glVertex3f(vertContainer_std_Lcontour_Leye.at(i)[0] - panel_width,					vertContainer_std_Lcontour_Leye.at(i)[1],			vertContainer_std_Lcontour_Leye.at(i)[2]); //1

		}	

		glVertex3f(vertContainer_std_Lcontour_Leye.at(n-1)[0],		vertContainer_std_Lcontour_Leye.at(n-1)[1] + panel_height_extra,			vertContainer_std_Lcontour_Leye.at(n-1)[2]); //0
		glVertex3f(vertContainer_std_Lcontour_Leye.at(n-1)[0] - panel_width,					vertContainer_std_Lcontour_Leye.at(n-1)[1] + panel_height_extra,			vertContainer_std_Lcontour_Leye.at(n-1)[2]); //1

		glEnd();
		glPopMatrix();
	}
	
}


void drawPanel_Reye(bool isStandard, double displayDist, double dispDepth){

	float panel_width = 40;
	float panel_height_extra = 20;
	int n = int(vertContainer_std_Rcontour_Reye.size());

	if(n > 0){
		glPushMatrix();

		glLoadIdentity();
		glTranslated(0, 0, displayDist - dispDepth + 2);
		glColor3f(0.0f, 0.0f, 0.0f);

		glBegin(GL_QUAD_STRIP);

		glVertex3f(vertContainer_std_Rcontour_Reye.at(0)[0] + panel_width,		vertContainer_std_Rcontour_Reye.at(0)[1] - panel_height_extra,			vertContainer_std_Rcontour_Reye.at(0)[2]); //0
		glVertex3f(vertContainer_std_Rcontour_Reye.at(0)[0],					vertContainer_std_Rcontour_Reye.at(0)[1] - panel_height_extra,			vertContainer_std_Rcontour_Reye.at(0)[2]); //1

		for (int i = 0; i < n; i++)
		{	
			glVertex3f(vertContainer_std_Rcontour_Reye.at(i)[0] + panel_width,		vertContainer_std_Rcontour_Reye.at(i)[1],			vertContainer_std_Rcontour_Reye.at(i)[2]); //0
			glVertex3f(vertContainer_std_Rcontour_Reye.at(i)[0],					vertContainer_std_Rcontour_Reye.at(i)[1],			vertContainer_std_Rcontour_Reye.at(i)[2]); //1

		}	

		glVertex3f(vertContainer_std_Rcontour_Reye.at(n-1)[0] + panel_width,		vertContainer_std_Rcontour_Reye.at(n-1)[1] + panel_height_extra,			vertContainer_std_Rcontour_Reye.at(n-1)[2]); //0
		glVertex3f(vertContainer_std_Rcontour_Reye.at(n-1)[0],					vertContainer_std_Rcontour_Reye.at(n-1)[1] + panel_height_extra,			vertContainer_std_Rcontour_Reye.at(n-1)[2]); //1

		glEnd();


		glBegin(GL_QUAD_STRIP);

		glVertex3f(vertContainer_std_Lcontour_Reye.at(0)[0],		vertContainer_std_Lcontour_Reye.at(0)[1] - panel_height_extra,			vertContainer_std_Lcontour_Reye.at(0)[2]); //0
		glVertex3f(vertContainer_std_Lcontour_Reye.at(0)[0] - panel_width,					vertContainer_std_Lcontour_Reye.at(0)[1] - panel_height_extra,			vertContainer_std_Lcontour_Reye.at(0)[2]); //1

		for (int i = 0; i < n; i++)
		{	
			glVertex3f(vertContainer_std_Lcontour_Reye.at(i)[0],		vertContainer_std_Lcontour_Reye.at(i)[1],			vertContainer_std_Lcontour_Reye.at(i)[2]); //0
			glVertex3f(vertContainer_std_Lcontour_Reye.at(i)[0] - panel_width,					vertContainer_std_Lcontour_Reye.at(i)[1],			vertContainer_std_Lcontour_Reye.at(i)[2]); //1

		}	

		glVertex3f(vertContainer_std_Lcontour_Reye.at(n-1)[0],		vertContainer_std_Lcontour_Reye.at(n-1)[1] + panel_height_extra,			vertContainer_std_Lcontour_Reye.at(n-1)[2]); //0
		glVertex3f(vertContainer_std_Lcontour_Reye.at(n-1)[0] - panel_width,					vertContainer_std_Lcontour_Reye.at(n-1)[1] + panel_height_extra,			vertContainer_std_Lcontour_Reye.at(n-1)[2]); //1

		glEnd();

		glPopMatrix();
	}
}

float adjustAmbient(double textDepth){
	float newAmbient = 0.4 - (textDepth - 20) * 0.005;
	return newAmbient;
}

void drawVertices(bool isStandard, int texNum, double displayDist, double dispDepth) {

	
    turnLightOn();

	// enable matrices for use in drawing below
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_BLEND);
	glEnable(GL_TEXTURE_2D);

	
	// bind the texture

	glBindTexture(GL_TEXTURE_2D, loaded_textures[texNum]);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	// activate and specify pointer to vertex array
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glPushMatrix();
	glLoadIdentity();
	glTranslated(0, 0, displayDist - dispDepth);

	if(isStandard){

		glVertexPointer(3, GL_FLOAT, 0, &vertices_vec_std[0]);
		glTexCoordPointer(2, GL_FLOAT, 0, &texcoors_vec_std[0]);
		glNormalPointer(GL_FLOAT, 0, &normals_vec_std[0]); //
		glColorPointer(3, GL_FLOAT, 0, &colors_vec_std[0]);
		glDrawElements(GL_TRIANGLES, indices_draw_triangle_vec_std.size(), GL_UNSIGNED_INT, &indices_draw_triangle_vec_std[0]);
	
	}else{

		glVertexPointer(3, GL_FLOAT, 0, &vertices_vec_cmp[0]);
		glTexCoordPointer(2, GL_FLOAT, 0, &texcoors_vec_cmp[0]);
		glNormalPointer(GL_FLOAT, 0, &normals_vec_cmp[0]); //
		glColorPointer(3, GL_FLOAT, 0, &colors_vec_cmp[0]);
		glDrawElements(GL_TRIANGLES, indices_draw_triangle_vec_cmp.size(), GL_UNSIGNED_INT, &indices_draw_triangle_vec_cmp[0]);

	}

	glPopMatrix();

	// deactivate vertex arrays after drawing
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

	glDisable(GL_LIGHTING); 


	switch(panel_state){
	case no_aperture:
		break;

	case black_aperture:

		glColor3f(0.0f, 0.0f, 0.0f);

		drawPanels(isStandard, displayDist, dispDepth);


		break;

	case red_aperture:

		glColor3f(0.4f, 0.0f, 0.0f);
		
		drawPanels(isStandard, displayDist, dispDepth);

		break;


	}
	

}

void drawBlockingPanels(double pandelSeparation){

	double panel_w = 20, panel_h = 80;

	// left panel
	glDisable(GL_TEXTURE_2D);

	glBegin(GL_QUADS);
	glVertex3f(-pandelSeparation / 2 - panel_w,   panel_h / 2,  0.0f);
	glVertex3f(-pandelSeparation / 2,             panel_h / 2,  0.0f);
	glVertex3f(-pandelSeparation / 2,            -panel_h / 2,  0.0f);
	glVertex3f(-pandelSeparation / 2 - panel_w,  -panel_h / 2,  0.0f);
	glEnd();

	// right panel
	glBegin(GL_QUADS);
	glVertex3f(pandelSeparation / 2,              panel_h / 2,  0.0f);
	glVertex3f(pandelSeparation / 2 + panel_w,    panel_h / 2,  0.0f);
	glVertex3f(pandelSeparation / 2 + panel_w,   -panel_h / 2,  0.0f);
	glVertex3f(pandelSeparation / 2,             -panel_h / 2,  0.0f);
	glEnd();
}

void drawFixation(double displayDist) {
	// draws a small fixation cross at the center of the display
	glDisable(GL_TEXTURE_2D);
	glColor3f(1.0f, 0.0f, 0.0f);
	glLineWidth(2.f);

	glPushMatrix();
	glLoadIdentity();
	glTranslated(0, 0, displayDist);
	double cross_length = 5;
	glBegin(GL_LINES);
	glVertex3d(cross_length / 2, 0, 0);
	glVertex3d(-cross_length / 2, 0, 0);
	glVertex3d(0, -cross_length / 2. , 0);
	glVertex3d(0, cross_length / 2. , 0);
	glEnd();
	glPopMatrix();

}

// This function seems to be used to shut down the system after use
void shutdown(){
	cout << "shutting down" << endl;
	responseFile.close(); // close this object
	if(resetScreen_betweenRuns){
		homeEverything(5000,4500);}
	cleanup();
	exit(0);
}
void cleanup() 
{
// Stop the optotrak
    optotrak->stopCollection();
    delete optotrak;
}
void initProjectionScreen(double _focalDist, const Affine3d &_transformation, bool synchronous)
{
	focalDistance = _focalDist;	
    screen.setWidthHeight(SCREEN_WIDE_SIZE, SCREEN_WIDE_SIZE*SCREEN_HEIGHT/SCREEN_WIDTH);
    screen.setOffset(alignmentX,alignmentY);
    screen.setFocalDistance(_focalDist);
    screen.transform(_transformation);
    cam.init(screen);
	if ( synchronous )
		moveScreenAbsolute(_focalDist,homeFocalDistance,4500);
	else
		moveScreenAbsoluteAsynchronous(_focalDist,homeFocalDistance,4500);
}
// Initialize Optotrak for use in the experiment
void initOptotrak()
{
    optotrak=new Optotrak2(); //intiailize the Optotrak object
    optotrak->setTranslation(calibration);

	//define Optotrak-specific variables
    int numMarkers=22;
    float frameRate=85.0f;
    float markerFreq=4600.0f;
    float dutyCycle=0.4f;
    float voltage = 7.0f;

	// run the intiailization method for the Optotrak, checking to see if ever (if == 0) and catch the error if so
    if ( optotrak->init("C:/cncsvisiondata/camerafiles/Aligned20111014",numMarkers, frameRate, markerFreq, dutyCycle,voltage) != 0)
    {   cerr << "Something during Optotrak initialization failed, press ENTER to continue. A error log has been generated, look \"opto.err\" in this folder" << endl;
        cin.ignore(1E6,'\n');
        exit(0);
    }

    // Read 10 frames of coordinates and fill the markers vector
    for (int i=0; i<10; i++)
    {
        updateTheMarkers();
    }
}
// run a method to define a vector that holds marker positions 
void updateTheMarkers()
{
	optotrak->updateMarkers();
	markers = optotrak->getAllMarkers();

}
// Initialize motors for moving screen around
void initMotors()
{
	//specify the speed for (objects,screen)
		if(resetScreen_betweenRuns){
		homeEverything(5000,4500);}

}

// Method that initializes the openGL parameters needed for creating the stimuli. 
// seems like this is not changed for each experiment (maybe for different experimental setup eg monitor)
void initRendering()
{   

	glClearColor(0.0,0.0,0.0,1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	/* Set depth buffer clear value */
	glClearDepth(1.0);
	/* Enable depth test */
	glEnable(GL_DEPTH_TEST);
	/* Set depth function */
	glDepthFunc(GL_LEQUAL);
	// scommenta solo se vuoi attivare lo shading degli stimoli

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glLineWidth(1.5);

	//texture-only
	glEnable(GL_MULTISAMPLE);
	glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);


}


void initVariables()
{
	interoculardistance = atof(parameters.find("IOD").c_str());
	display_distance = str2num<double>(parameters.find("dispDepth"));

	// eye coordinates

	eyeRight = Vector3d(interoculardistance / 2, 0, 0);
	eyeLeft = Vector3d(-interoculardistance / 2, 0, 0);

	eyeMiddle = Vector3d(0, 0, 0);

	stimulus_height = tan((DEG2RAD * visual_angle)/2) * 2 * (abs(display_distance));

	stimulus_bgwidth = ratio_bgwidth_height * stimulus_height;
	stimulus_visiblewidth = ratio_visiblewidth_height * stimulus_height;
	stimulus_width = stimulus_bgwidth;

	//turnLightOn();


}

void turnLightOn(){

	glShadeModel(GL_SMOOTH); // enable Smooth Shading
	glEnable(GL_LIGHTING); // enable lighting
	glEnable(GL_LIGHT1); 
	glEnable(GL_NORMALIZE); //so we don't need to normalize our normal for surfaces

	glPushMatrix();
	glLoadIdentity();

	// Light source parameters
	GLfloat LightAmbient[]= { amb_intensity, 0.0f, 0.0f, 1.0f }; // non-directional & overall light (r,g,b,alpha): dark part
	GLfloat LightDiffuse[]= { max_intensity - amb_intensity, 0.0f, 0.0f, 1.0f }; // light created by the light source (directional light; r,g,b,alpha): bright part
	GLfloat LightPosition[]= { 0.0f, 1.f, lightDir_z, 0.0f }; // Light Position (x, y, z, 1.0f); if w==0, directional; if w==1, positional lights. Attenuation can be applied only to the positional light 
	
	//setting the light
	glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient); //setup the ambient light
	glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse); //setup the diffuse light
	glLightfv(GL_LIGHT1, GL_POSITION,LightPosition); //position the light

	glPopMatrix();

	//Attenuation: based on the distance between the light source & the object surface
	//glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, 1.0); //default 1.0
	//glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, atten_weight*2); //default 0.0
	//glLightf(GL_LIGHT1, GL_QUADRATIC_ATTENUATION, atten_weight); //default 0.0

}

void initBlock()
{
	// initialize the trial matrix
	trial.init(parameters);
	//trial.next();
	trial.next(false);

	trialNum = 1;
}

// Initialize the streams, open the file and write to it
void initStreams()
{
	// Initialize the parameter file starting from the file parameters.txt, if the file does not exist, it tells you
	
	//parametersFileName = experiment_directory + "/parameters_summer22-ailin-hapticRemap-matching-staircase.txt";
	
	parametersFile.open(parametersFileName.c_str());
	parameters.loadParameterFile(parametersFile);


	// Subject name
	subjectName = parameters.find("SubjectName");

	string responseFileName = experiment_directory +"/"+ subjectName  + ".txt";
	// Principal streams files
	if (util::fileExists(experiment_directory +"/"+subjectName+".txt") && subjectName != "junk")
	{
		string error_on_file_io = experiment_directory+"/"+subjectName+".txt" + string(" already exists");
		cerr << error_on_file_io << endl;
		MessageBox(NULL, (LPCSTR)"FILE ALREADY EXISTS\n Please check the parameters file.",NULL, NULL);
		shutdown();
	}
	responseFile.open(responseFileName.c_str());
	responseFile << fixed << responseFile_headers << endl;
	
	//globalTimer.start();
}


//Central function for projecting image onto the screen
void drawGLScene()
{
		glDrawBuffer(GL_BACK_RIGHT);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glClearColor(0.0,0.0,0.0,1.0);
		cam.setEye(eyeLeft);

		drawStimulus();
		//drawPanel_Leye(true, display_distance_jittered_std, depth_std_disp);
		drawInfo();
		
		
		// Draw right eye view
		glDrawBuffer(GL_BACK_LEFT);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glClearColor(0.0,0.0,0.0,1.0);
		cam.setEye(eyeRight); 
		

		drawStimulus();
		//drawPanel_Reye(true, display_distance_jittered_std, depth_std_disp);
		drawInfo();
		

	glutSwapBuffers();
}

void update(int value)
{
    glutPostRedisplay();
    glutTimerFunc(TIMER_MS, update, 0);
}

void drawProgressBar() {
	
	glPushMatrix();
	glLoadIdentity();
	glTranslated(0, 0, display_distance);

	glColor3f(0.2, 0.2, 0.2);
	glBegin(GL_LINE_LOOP);
	glVertex3f(-50, 5, 0);
	glVertex3f(50, 5, 0);
	glVertex3f(50, -5, 0);
	glVertex3f(-50, -5, 0);
	glEnd();

	glColor3f(0.1, 0.3, 0.1);
	glBegin(GL_POLYGON);
	glVertex3f(-50, 5, 0);
	glVertex3f(-50 + percentComplete, 5, 0);
	glVertex3f(-50 + percentComplete, -5, 0);
	glVertex3f(-50, -5, 0);
	glEnd();
	glPopMatrix();
}


void drawInfo()
{
	// displays relevant information to the screen
	if ( visibleInfo )
	{
		glDisable(GL_COLOR_MATERIAL);
		glDisable(GL_BLEND);
		GLText text;
		if (gameMode)
			text.init(SCREEN_WIDTH, SCREEN_HEIGHT, glWhite, GLUT_BITMAP_HELVETICA_18);
		else
			text.init(640, 480, glWhite, GLUT_BITMAP_HELVETICA_12);

	
		text.enterTextInputMode();

		switch(current_stage){
			case stimulus_preview:
				glColor3fv(glRed);
				text.draw("Welcome! press + to start training");
				text.draw("# Name: " + subjectName);
				text.draw("# IOD: " + stringify<double>(interoculardistance));
				text.draw("# STD depth texture: " + stringify<double>(depth_std_text));
				text.draw("# STD depth stereo: " + stringify<double>(depth_std_disp));
				text.draw("                           ");
				text.draw("# LIGHT amb intensity: " + stringify<double>(amb_intensity));
				text.draw("# LIGHT z: " + stringify<float>(lightDir_z));
				
				
				//text.draw("# current enum id: " + stringify<int>(panel_state));
				
				
					// check if mirror is calibrated
				if (abs(mirrorAlignment - 45.0) > 0.2)
					text.draw("# !!!!Mirror Alignment = " + stringify<double>(mirrorAlignment));

				break;

			case trial_respond:

				glColor3fv(glRed);
				text.draw("# Name: " + subjectName);
				text.draw("# IOD: " + stringify<double>(interoculardistance));
				text.draw("# trial: " + stringify<int>(trialNum));
				//text.draw("# time: " + stringify<double>(ElapsedTime));
				text.draw("# STD depth texture: " + stringify<double>(depth_std_text));
				text.draw("# STD depth stereo: " + stringify<double>(depth_std_disp));
				text.draw("                           ");
				text.draw("# CMP depth: " + stringify<double>(depth_cmp));

									// check if mirror is calibrated
				if (abs(mirrorAlignment - 45.0) > 0.2)
					text.draw("# !!!!Mirror Alignment = " + stringify<double>(mirrorAlignment));
				break;

			case break_time:
				if(training && trialNum >= trainNum_cap){
					text.draw("Pls call the experimenter!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
				}else{
					glColor3fv(glRed);
					text.draw("Break time! Press + to continue");
				}
				break;

			case exp_completed:
				glColor3fv(glWhite);
				text.draw("The experiment is over. Thank you! :)");
				break;
		}
		text.leaveTextInputMode();
		glEnable(GL_COLOR_MATERIAL);
		glEnable(GL_BLEND);
	}
}


// Funzione che gestisce il ridimensionamento della finestra
void handleResize(int w, int h)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glViewport(0,0,SCREEN_WIDTH, SCREEN_HEIGHT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

}

void drawStimulus()
{
	//enum Stages {stimulus_preview, prep_trial, trial_fixate_first, trial_present_first,
	//trial_fixate_second, trial_present_second, trial_respond, break_time, exp_completed};
	switch(current_stage){

		case stimulus_preview:
			// testing texture, build vertices for curved surface
			drawVertices(true, texnum_std, display_distance_jittered_std, depth_std_disp);
			
			break;

		case trial_fixate_first:
			if(std_is_first){
				drawFixation(display_distance_jittered_std);
			}else{
				drawFixation(display_distance_jittered_cmp);
			}
			break;

		case trial_present_first:
			if(std_is_first)
				drawVertices(true, texnum_std, display_distance_jittered_std, depth_std_disp);
			else
				drawVertices(false, texnum_cmp, display_distance_jittered_cmp, depth_cmp);

			break;

		case trial_fixate_second:
			if(std_is_first){
				drawFixation(display_distance_jittered_cmp);
			}else{
				drawFixation(display_distance_jittered_std);
			}
			break;

		case trial_present_second:
			if(std_is_first)
				drawVertices(false, texnum_cmp, display_distance_jittered_cmp, depth_cmp);
			else
				drawVertices(true, texnum_std, display_distance_jittered_std, depth_std_disp);
			break;

		case break_time:
			drawProgressBar();
			break;

	}
}




void initTrial()
{
	current_stage = prep_trial;
	initProjectionScreen(display_distance);
	trial_timer.reset();
	trial_timer.start();

	if(rand() % 2){
		std_is_first = true;
	}else{
		std_is_first = false;
	}

	if(training){

		depth_std = rand() % depth_training_range + depth_training_min;
		depth_std_disp = depth_std;
		depth_std_text = depth_std;
		depth_cmp = rand() % depth_training_range + depth_training_min;
	}
	else {

		stairID = trial.getCurrent().second->getCurrentStaircase()->getID();
		stair_reversal = trial.getCurrent().second->getCurrentStaircase()->getReversals();
		ascending = trial.getCurrent().second->getCurrentStaircase()->getAscending();

		depth_std = trial.getCurrent().first["stdDepth"];
		depth_stdDelta = trial.getCurrent().first["stdDepthDelta"];
		depth_std_text = depth_std + depth_stdDelta;
		depth_std_disp = depth_std - depth_stdDelta;

		depth_cmp = trial.getCurrent().second->getCurrentStaircase()->getState();
	}
	
	jitter_z_std = (rand() % 201 - 100)/10.0; // from -10 to 10
	display_distance_jittered_std = display_distance + jitter_z_std;

	jitter_z_cmp = (rand() % 201 - 100)/10.0; // from -10 to 10
	display_distance_jittered_cmp = display_distance + jitter_z_cmp;

	normalizer_to_uv_std = normalizer_to_uv_base + ((rand() % 13) - 6); // from base-6 to base+6
	normalizer_to_uv_cmp = normalizer_to_uv_base + ((rand() % 13) - 6); // from base-6 to base+6

	if(abs(depth_std_text - depth_std_disp) < 0.1){
		buildVertices_congruent(true, depth_std_text, normalizer_to_uv_std);
	}else{
		buildVertices_incongruent(true, depth_std_text, depth_std_disp, display_distance_jittered_std, normalizer_to_uv_std);
	}		
	buildVertices_congruent(false, depth_cmp, normalizer_to_uv_cmp);


	texnum_std = rand() % 50 + 1;
	texnum_cmp = rand() % 50 + 1;
	if(texnum_cmp == texnum_cmp){
		texnum_cmp = rand() % 50 + 1;
	}

	trial_timer.reset();				
	trial_timer.start();
	ElapsedTime = 0;

	current_stage = trial_fixate_first;	
}

void onlineTrial(){
	//stimulus_preview, prep_trial, trial_fixate_first, trial_present_first, trial_fixate_second, trial_present_second, trial_respond,

	switch(current_stage){

	case trial_fixate_first:

		if (ElapsedTime > fixateTime) {
			current_stage = trial_present_first;
		}
		break;

	case trial_present_first:

		if (ElapsedTime > (fixateTime + viewTime)) {
			current_stage = trial_fixate_second;
		}
		break;

	case trial_fixate_second:

		if (ElapsedTime > (2 * fixateTime + viewTime)) {
			current_stage = trial_present_second;
		}
		break;

	case trial_present_second:

		if (ElapsedTime > (2 * fixateTime + 2 * viewTime) ){
			current_stage = trial_respond;
		}
		break;
	}

}

void advanceTrial()
{
	//subjName\tIOD\tblockN\ttrialN\tdisplayDistance\tvisualAngle\tclyHorizontal\ttexnum\ttextNomralizer\ttestDepth\tprobeStart\tprobeDepth\ttime
	if(training){
		if(trialNum < trainNum_cap){
			beepOk(4);
			trialNum++;
			initTrial();
		}else{
			beepOk(2);
			trialNum++;
			current_stage = break_time;
			visibleInfo = true;			
		}

	}else{
		//subjName\tIOD\tblockN\ttrialN\tdisplayDistance\tvisualAngle\tshapeID\ttexnum\ttextNomralizer\ttestDepth\tprobeDepthInit\tprobeDepth\tRT
		responseFile << fixed <<
		subjectName << "\t" <<		
		interoculardistance << "\t" <<
		blkNum << "\t" <<
		trialNum << "\t" <<
		display_distance << "\t" <<
		visual_angle << "\t" <<
		stimulus_height << "\t" <<
		texnum_std << "\t" <<
		normalizer_to_uv_std << "\t" <<
		texnum_cmp << "\t" <<
		normalizer_to_uv_cmp << "\t" <<
		depth_std_text << "\t" <<
		depth_std_disp << "\t" <<
		depth_cmp << "\t" <<
		std_is_first << "\t" <<
		respond_first_deeper << "\t" <<
		respond_cmp_deeper << "\t" <<
		ElapsedTime  << "\t" <<
		stairID << "\t" <<
		stair_reversal << "\t" <<
		ascending << endl;	

		trial.next(respond_cmp_deeper);	

		if (!trial.isEmpty()){

			if(trialNum % 24 == 0){
				beepOk(4);				
				percentComplete = trialNum /(totalTrNum/100.0);			
				current_stage = break_time;
				visibleInfo = true;

			}else{
				beepOk(4);
				trialNum++;
				initTrial();
			}
			
		}else{
			beepOk(1);
			responseFile.close();
			visibleInfo=true;
			current_stage = exp_completed;
		}

	}


}



void handleKeypress(unsigned char key, int x, int y)
{   
	switch (key){ // key presses that work regardless of the stage
	
		case 27:	//corrisponde al tasto ESC
			shutdown();		
		break;

		case 'i':
			visibleInfo=!visibleInfo;
		break;

		case 'o':
			panel_state = panelStates((panel_state + 1)%3);
			break;

		case '1':
			if(current_stage == stimulus_preview){
				if (depth_std_text > depth_inc)
					depth_std_text = depth_std_text - depth_inc;

				amb_intensity = adjustAmbient(depth_std_text);
				initPreviewStimulus(depth_std_text, depth_std_disp);
			}else if(current_stage == trial_respond){
				trial_timer.stop();

				respond_first_deeper = true;
				if(std_is_first){
					respond_cmp_deeper = false;
				}else{
					respond_cmp_deeper = true;
				}

				advanceTrial();
			}
			break;

		case '2':
			if(current_stage == stimulus_preview){

				depth_std_text = depth_std_text + depth_inc;
				amb_intensity = adjustAmbient(depth_std_text);
				initPreviewStimulus(depth_std_text, depth_std_disp);
			}else if(current_stage == trial_respond){
				trial_timer.stop();

				respond_first_deeper = false;
				if(std_is_first){
					respond_cmp_deeper = true;
				}else{
					respond_cmp_deeper = false;
				}

				advanceTrial();
			}
			break;




		case '+':
			switch(current_stage){
				case stimulus_preview:
					
					beepOk(5);
					visibleInfo = false;
					initBlock();
					initTrial();
					
				break;

				case trial_respond:
					trial_timer.stop();
					advanceTrial();
				break;

				case break_time:
					beepOk(5);
					visibleInfo = false;
					trialNum++;
					initTrial();
				break;		
			}
			break;

		case 'T':
		case 't':
			if(training){
				training = false;
				beepOk(6);
				trialNum = 0;
				current_stage = break_time;
			}
			break;

		///////////////////////////////////////////////
		// key presses for adjusting stimulus previe

		case '7':
			if(current_stage == stimulus_preview){
				
				if (depth_std_disp > depth_inc)
					depth_std_disp = depth_std_disp - depth_inc;

				initPreviewStimulus(depth_std_text, depth_std_disp);
				
			}
			break;

		case '8':
			if(current_stage == stimulus_preview){
				
				depth_std_disp = depth_std_disp + depth_inc;

				initPreviewStimulus(depth_std_text, depth_std_disp);
				
			}
			break;



			
	}

}

/***** SOUNDS *****/
void beepOk(int tone)
{

	switch(tone)
	{

	case 1: //high pitch beep
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-8_lowpass.wav", NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 4: //mellow and good for trials
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-440-pluck.wav", NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 2: //reject like buzz
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-10.wav", NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 3: //reject short
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-reject.wav", NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 5: //mellow and good for trials
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-highBubblePop.wav", NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 6: //mellow and good for trials
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-lowBubblePop.wav", NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 15:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-rising.wav",
			NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 16:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-falling.wav",
			NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 17:
	PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-440-pluck-5below.wav",
		NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 18: // light click
	PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-click3MS.wav",
		NULL, SND_FILENAME | SND_ASYNC);
	break;

	}
	return;
}

void idle()
{
	onlineTrial();
	online_apparatus_alignment();
	ElapsedTime = trial_timer.getElapsedTimeInMilliSec();
}

/*** Online operations ***/
void online_apparatus_alignment()
{
	updateTheMarkers();
    // mirror alignment check
    mirrorAlignment = asin(
        abs((markers.at(mirror1).p.z() - markers.at(mirror2).p.z())) /
        sqrt(
        pow(markers.at(mirror1).p.x() - markers.at(mirror2).p.x(), 2) +
        pow(markers.at(mirror1).p.z() - markers.at(mirror2).p.z(), 2)
        )
        ) * 180 / M_PI;

    // screen Y alignment check
    screenAlignmentY = asin(
		abs((markers.at(screen1).p.y() - markers.at(screen3).p.y())) /
		sqrt(
		pow(markers.at(screen1).p.x() - markers.at(screen3).p.x(), 2) +
		pow(markers.at(screen1).p.y() - markers.at(screen3).p.y(), 2)
		)
		) * 180 / M_PI;

    // screen Z alignment check
    screenAlignmentZ = asin(
		abs(markers.at(screen1).p.z() - markers.at(screen2).p.z()) /
		sqrt(
		pow(markers.at(screen1).p.x() - markers.at(screen2).p.x(), 2) +
		pow(markers.at(screen1).p.z() - markers.at(screen2).p.z(), 2)
		)
		) * 180 / M_PI *
		abs(markers.at(screen1).p.x() - markers.at(screen2).p.x()) /
		(markers.at(screen1).p.x() - markers.at(screen2).p.x());

}

int LoadGLTextures()  // Load PNG And Convert To Textures
{

	for (int i = 1; i <= 50; i++) {
		std::stringstream ss;
		ss << i;

		string texturePath = "fall16-ailin-stimulusTest/0/polkadots" + ss.str() + ".png";
		 loaded_textures[i] = SOIL_load_OGL_texture
		(
			texturePath.c_str(),
			SOIL_LOAD_AUTO,
			SOIL_CREATE_NEW_ID,
			SOIL_FLAG_MULTIPLY_ALPHA
		);
	}


	return true; // Return Success
}

// this is run at compilation because it's titled 'main'
int main(int argc, char*argv[])  
{
	//functions from cncsvision packages
	mathcommon::randomizeStart();
		
	// initializing glut (to use OpenGL)
	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_STEREO | GLUT_MULTISAMPLE);

	glutGameModeString("1024x768:32@85"); //resolution  
	glutEnterGameMode();
	glutFullScreen();
	
	// initializes optotrak and velmex motors
	initOptotrak();
	initMotors();
	
	initRendering(); // initializes the openGL parameters needed for creating the stimuli
	
	initStreams(); // streams as in files for writing data

	initVariables();

	LoadGLTextures();

	initPreviewStimulus(depth_std_text, depth_std_disp);

	initProjectionScreen(display_distance);

	// glut callback, OpenGL functions that are infinite loops to constantly run 

	glutDisplayFunc(drawGLScene); // keep drawing the stimuli

	glutKeyboardFunc(handleKeypress); // check for keypress

	glutReshapeFunc(handleResize);

	glutIdleFunc(idle);

	glutTimerFunc(TIMER_MS, update, 0);

	glutSetCursor(GLUT_CURSOR_NONE);

	//boost::thread initVariablesThread(&initVariables); 

	glutMainLoop();

	//cleanup();
	return 0;
}
