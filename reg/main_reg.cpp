////gdal
#include "gdal_priv.h"
#include "gdalwarper.h"
#include "ogr_core.h"
#include "ogr_api.h"
#include "cpl_conv.h" // for CPLMalloc()
#include <ogrsf_frmts.h>

#include "iostream"
#include <algorithm>
#include <sstream>
#include <iterator>
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#include <experimental/filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <ICP.h>
#include <strutil.h>
#include <math.h>
#include <Eigen/Geometry>
#include <numeric>
#include <argparse.hpp>
#include <nanoflann.hpp>
#include <chrono>
using namespace std;
namespace fs = std::experimental::filesystem;
using namespace std::chrono;
#define PI 3.14159265;

///////////////////////////////////////////////////////////////////////////////

template <typename T>
struct PointCloud
{
    struct Point
    {
        T  x, y, z;
    };

    std::vector<Point>  pts;

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return pts.size(); }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate value, the
    //  "if/else's" are actually solved at compile time.
    inline T kdtree_get_pt(const size_t idx, const size_t dim) const
    {
        if (dim == 0) return pts[idx].x;
        else if (dim == 1) return pts[idx].y;
        else return pts[idx].z;
    }
    // Optional bounding-box computation: return false to default to a standard bbox computation loop.
//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};

class TILE_INFO {
public:
    TILE_INFO() {
        this->utm_bbox[0] = 0;
        this->utm_bbox[1] = 0;
        this->utm_bbox[2] = 0;
        this->utm_bbox[3] = 0;
        this->ref_var = 0;
        this->src_var = 0;
        this->tile_id = 0;
        this->src_lambda1 = 0;
        this->src_lambda2 = 0;
        this->src_R=0;
        this->ref_lambda1 = 0;
        this->ref_lambda2 = 0;
        this->ref_R = 0;
    }
    double utm_bbox[4];
    double ref_var;
    double src_var;
    double ref_valid;
    double src_valid;
    double src_lambda1;
    double src_lambda2;
    double src_R;
    double ref_lambda1;
    double ref_lambda2;
    double ref_R;
    int tile_id;
    double rotation[9];
    double translation[3];
};

class CORR {
public:
    CORR() {};
    CORR(gs::Point src,gs::Point ref,double weight) {
        this->src = src;
        this->ref = ref;
        this->w = weight;
    }
    CORR(gs::Point src, gs::Point ref,std::vector<float> src_norm, std::vector<float> ref_norm, double weight) {
        this->src = src;
        this->ref = ref;
        this->w = weight;
        this->src_norm = src_norm;
        this->ref_norm = ref_norm;
    }
    gs::Point src;
    gs::Point ref;
    vector<float> src_norm;
    vector<float> ref_norm;
    double w;
};

struct Parameters {
    std::string type_ = "rigid";
    std::string ref_path_, src_path_, out_dir_;
    double valid_ratio_threshold_ = 0.5;
    double var_min_=4.0;
    int verbose_=1;
    // for General 
    double outlier_thresh_ = -1;
    double src_resolution_x_ = 0.5;
    double src_resolution_y_ = 0.5;
    double ref_resolution_x_ = 0.5;
    double ref_resolution_y_ = 0.5;
    bool plane_or_not_=false;
    bool multi_or_not_=false;
    bool brute_force_search_or_not_=false;
    double xmin_m_=-3; // for brute-force search bounding
    double xmax_m_=3;
    double ymin_m_=-3;
    double ymax_m_=3;
    double rough_step_m_=0.5;
    int rough_step_min_resolution_=3; //local_minimal_within_range=rough_step_min_resolution_*rough_step_m_
    int rough_num_mins_=3; // number of local minimal
    double fine_step_m_=0.5;
    int gen_data_ = 0;

    // for rough align
    double icp_num_pts_ratio_ = 0.01;
    int rough_num_pts_ = 24000;
    int search_half_size_pixel_ = 5;
    double valid_ratio_threshold_rough_=0.5;
    double rough_icp_max_iter_=200;
    double rough_icp_rmse_threshold_=1e-5;
    double trimmed_ratio_=0.9;
    double rough_icp_percentage_threshold_=1e-5;

    // for fine align
    double tilesize_m_ = 100;
    double HARRIS_step_=3;
    double HARRIS_k_=0.05;
    int n_blocks_col_=5;
    int n_blocks_row_=5;
    int num_filter_tiles_=100;
    double fine_icp_max_iter_=200;
    double fine_icp_rmse_threshold_=1e-5;
    double fine_icp_percentage_threshold_=1e-5;
    double fine_search_half_pixels_=10;
    int fine_resolution_factor_=5;
};

void parse_shp(string shp_path, vector<string>& tif_paths, vector<double>& tif_bboxs) {
    GDALDataset* shps;
    shps = (GDALDataset*)GDALOpenEx(shp_path.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
    int a = shps->GetLayerCount();
    OGRLayer* l1 = shps->GetLayer(0);
    OGRFeature* poFeature;
    l1->ResetReading();
    while ((poFeature = l1->GetNextFeature()) != NULL)
    {
        int c = poFeature->GetFieldCount();
        OGRGeometry* shape = poFeature->GetGeometryRef();
        OGREnvelope psEnvelope;
        shape->getBoundary()->getEnvelope(&psEnvelope);
        auto&& field = (*poFeature)[0];
        string tif_path = field.GetAsString();
        tif_bboxs.push_back(psEnvelope.MinX);
        tif_bboxs.push_back(psEnvelope.MaxX);
        tif_bboxs.push_back(psEnvelope.MinY);
        tif_bboxs.push_back(psEnvelope.MaxY);
        tif_paths.push_back(tif_path);
    }
    GDALClose(shps);
}

void export_tfw(string path, double* transform) {
    std::ofstream out(path.c_str());
    out << std::fixed << std::setprecision(8) << transform[1] << endl;
    out << std::fixed << std::setprecision(8) << transform[2] << endl;
    out << std::fixed << std::setprecision(8) << transform[4] << endl;
    out << std::fixed << std::setprecision(8) << transform[5] << endl;
    out << std::fixed << std::setprecision(8) << transform[0] << endl;
    out << std::fixed << std::setprecision(8) << transform[3] << endl;
    out.close();

}

void get_aoi_bbox_pixels(double* aoi_utm_bbox, double* ref_bbox, double* src_bbox, int ref_width, int ref_height, int src_width, int src_height,double ref_reso,double src_reso,
    int* aoi_ref_bbox, int* aoi_src_bbox) {
    aoi_ref_bbox[0] = floor((aoi_utm_bbox[0] - ref_bbox[0]+ref_reso/2) / ref_reso);
    aoi_ref_bbox[1] = floor((aoi_utm_bbox[1] - ref_bbox[0]+ ref_reso/2) / ref_reso);
    aoi_ref_bbox[2] = floor((ref_bbox[3]-aoi_utm_bbox[2]+ ref_reso/2) / ref_reso);
    aoi_ref_bbox[3] = floor((ref_bbox[3]-aoi_utm_bbox[3]+ ref_reso/2) / ref_reso);
    aoi_src_bbox[0] = floor((aoi_utm_bbox[0] - src_bbox[0]+ src_reso/2) /src_reso);
    aoi_src_bbox[1] = floor((aoi_utm_bbox[1] - src_bbox[0]+ src_reso/2) / src_reso);
    aoi_src_bbox[2] = floor((src_bbox[3]-aoi_utm_bbox[2]+ src_reso/2) / src_reso);
    aoi_src_bbox[3] = floor((src_bbox[3]-aoi_utm_bbox[3]+ src_reso/2) / src_reso);
    aoi_ref_bbox[0] = max(aoi_ref_bbox[0], 0);
    aoi_ref_bbox[0] = min(aoi_ref_bbox[0], ref_width - 1);
    aoi_ref_bbox[1] = max(aoi_ref_bbox[1], 0);
    aoi_ref_bbox[1] = min(aoi_ref_bbox[1], ref_width - 1);
    aoi_ref_bbox[2] = min(aoi_ref_bbox[2], ref_height - 1);
    aoi_ref_bbox[2] = max(aoi_ref_bbox[2], 0);
    aoi_ref_bbox[3] = min(aoi_ref_bbox[3], ref_height - 1);
    aoi_ref_bbox[3] = max(aoi_ref_bbox[3], 0);
    aoi_src_bbox[0] = max(aoi_src_bbox[0], 0);
    aoi_src_bbox[0] = min(aoi_src_bbox[0], src_width - 1);
    aoi_src_bbox[1] = max(aoi_src_bbox[1], 0);
    aoi_src_bbox[1] = min(aoi_src_bbox[1], src_width - 1);
    aoi_src_bbox[2] = min(aoi_src_bbox[2], src_height - 1);
    aoi_src_bbox[2] = max(aoi_src_bbox[2], 0);
    aoi_src_bbox[3] = min(aoi_src_bbox[3], src_height - 1);
    aoi_src_bbox[3] = max(aoi_src_bbox[3], 0);
}

void compute_tile_var(double* data, int width, int height,double& var, double& valid_ratio) {
    double mean=0;
    var = 0;
    int valid_cnt = 0;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            if (isnan(data[i * width + j])) {
                continue;
            }
            mean += data[i * width + j];
            valid_cnt++;
        }
    }
    mean /=valid_cnt;
    valid_ratio = (double)valid_cnt / (double)(width * height);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            if (isnan(data[i * width + j])) {
                continue;
            }
            var += pow(data[i * width + j]-mean,2.0);
        }
    }
    var/= valid_cnt;
    var = std::sqrt(var);
}

void write_ptcloud(double* data, int width, int height, double resolution, string out_path,double offx,double offy) {
    std::ofstream outfile(out_path);
    for (int row = 0; row < height; ++row) {
        for (int col = 0; col < width; ++col) {
            if (isnan(data[col + row * width])) {
                continue;
            }
            outfile << std::fixed << std::setw(11) << std::setprecision(11) << (double)col * resolution+offx << " " <<
                std::fixed << std::setw(11) << std::setprecision(11) << (double)(height-row) * resolution+offy << " " <<
                std::fixed << std::setw(11) << std::setprecision(11) << data[col + row * width] << "\n";
        }
    }
    outfile.close();
}

void write_ptcloud_corrs(std::vector<pair<gs::Point*, gs::Point*>>& corrs, string out_ref_path, string out_src_path) {
    std::ofstream out_ref(out_ref_path);
    std::ofstream out_src(out_src_path);
    for (int i = 0; i < corrs.size(); ++i) {
        out_src << std::fixed << std::setw(11) << std::setprecision(11) << (double)corrs[i].first->pos[0] << " " <<
            std::fixed << std::setw(11) << std::setprecision(11) << (double)corrs[i].first->pos[1] << " " <<
            std::fixed << std::setw(11) << std::setprecision(11) << (double)corrs[i].first->pos[2] << "\n";
    }
    for (int i = 0; i < corrs.size(); ++i) {
        out_ref << std::fixed << std::setw(11) << std::setprecision(11) << (double)corrs[i].second->pos[0] << " " <<
            std::fixed << std::setw(11) << std::setprecision(11) << (double)corrs[i].second->pos[1] << " " <<
            std::fixed << std::setw(11) << std::setprecision(11) << (double)corrs[i].second->pos[2] << "\n";
    }

    out_ref.close();
    out_src.close();
}

void write_ptcloud_gs(std::vector<gs::Point>& pts, string out_path,double offx,double offy) {
    std::ofstream outfile(out_path);
    for (int i= 0; i < pts.size(); ++i) {
        outfile << std::fixed << std::setw(11) << std::setprecision(11) << (double)pts[i].pos[0]+offx << " " <<
            std::fixed << std::setw(11) << std::setprecision(11) << (double)pts[i].pos[1]+offy << " " <<
            std::fixed << std::setw(11) << std::setprecision(11) << (double)pts[i].pos[2] << "\n";

    }
    outfile.close();
}

void convert_ptcloud(double* data, int width, int height,double reso_x,double reso_y,int resolution_factor, vector<gs::Point>& ptcloud, double offx, double offy) {
    ptcloud.clear();
    for (int row = 0; row < height; row+=resolution_factor) {
        for (int col = 0; col < width; col+=resolution_factor) {
            if (isnan(data[col + row * width])) {
                continue;
            }
            gs::Point pt((double)(col+0.5) * reso_x+offx, (double)(- row-0.5) * reso_y+offy, data[col + row * width]);
            ptcloud.push_back(pt);
        }
    }

}

void update_ptcloud(vector<gs::Point>& ptcloud, Eigen::Matrix4d T,double offx_in,double offy_in,double offx_rot,double offy_rot) {
    double tmp[3];
    for (int i = 0; i < ptcloud.size(); ++i) {
        ptcloud[i].pos[0] = ptcloud[i].pos[0] + offx_in - offx_rot;
        ptcloud[i].pos[1] = ptcloud[i].pos[1] + offy_in - offy_rot;
        Eigen::Vector4d pt(ptcloud[i].pos[0], ptcloud[i].pos[1], ptcloud[i].pos[2], 1);
        pt = T * pt;

        //gs::rotate_mat(ptcloud[i].pos, rot, tmp);
        ptcloud[i].pos[0] = pt(0);
        ptcloud[i].pos[1] = pt(1);
        ptcloud[i].pos[2] = pt(2);
        
        ptcloud[i].pos[0] = ptcloud[i].pos[0]+offx_rot - offx_in;
        ptcloud[i].pos[1] = ptcloud[i].pos[1]+offy_rot - offy_in;
    }
}

void update_ptcloud_offset(vector<gs::Point>& ptcloud, double offx_in, double offy_in, double offx_out, double offy_out) {
    double tmp[3];
    for (int i = 0; i < ptcloud.size(); ++i) {
        ptcloud[i].pos[0] = ptcloud[i].pos[0] + offx_in - offx_out;
        ptcloud[i].pos[1] = ptcloud[i].pos[1] + offy_in - offy_out;
    }
}

void load_rmse(std::string path, double** rmses) {
    ifstream ifs(path);
    int x, y;
    double rmse;
    char line[1000];
    string line_str;
    while (!ifs.eof()) {
        //ifs >> line;
        getline(ifs, line_str);
        sscanf(line_str.c_str(), "%d %d %lf", &x, &y, &rmse);
        rmses[y][x] = rmse;
    }
}

void write_shp(std::vector<TILE_INFO> tile_infos, string out_path) {
    GDALAllRegister();
    const char* pszDriverName = "ESRI Shapefile";
    GDALDriver* poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName(pszDriverName);
    if (poDriver == NULL)
    {
        printf("%s driver not available.\n", pszDriverName);
        exit(1);
    }
    GDALDataset* poDS;
    poDS = poDriver->Create(out_path.c_str(), 0, 0, 0, GDT_Unknown, NULL);
    if (poDS == NULL)
    {
        printf("Creation of output file failed.\n");
        exit(1);
    }
    OGRLayer* poLayer;

    poLayer = poDS->CreateLayer("tiles", NULL, wkbPolygon, NULL);
    if (poLayer == NULL)
    {
        printf("Layer creation failed.\n");
        exit(1);
    }
    OGRFieldDefn oField1("Ref_Var", OFTString);
    OGRFieldDefn oField2("Src_Var", OFTString);
    OGRFieldDefn oField3("Tile_ID", OFTString);
    OGRFieldDefn oField4("Ref_Valid", OFTString);
    OGRFieldDefn oField5("Src_Valid", OFTString);
    OGRFieldDefn oField6("Src_lmd1", OFTString);
    OGRFieldDefn oField7("Src_lmd2", OFTString);
    OGRFieldDefn oField8("Src_R", OFTString);
    OGRFieldDefn oField9("Ref_lmd1", OFTString);
    OGRFieldDefn oField10("Ref_lmd2", OFTString);
    OGRFieldDefn oField11("Ref_R", OFTString);
    oField1.SetWidth(32);
    oField2.SetWidth(32);
    oField3.SetWidth(32);
    oField4.SetWidth(32);
    oField5.SetWidth(32);
    oField6.SetWidth(32);
    oField7.SetWidth(32);
    oField8.SetWidth(32);
    oField9.SetWidth(32);
    oField10.SetWidth(32);
    oField11.SetWidth(32);

    if (poLayer->CreateField(&oField1) != OGRERR_NONE || poLayer->CreateField(&oField2) != OGRERR_NONE || 
        poLayer->CreateField(&oField3) != OGRERR_NONE || poLayer->CreateField(&oField4) != OGRERR_NONE || 
        poLayer->CreateField(&oField5) != OGRERR_NONE || poLayer->CreateField(&oField6) != OGRERR_NONE ||
        poLayer->CreateField(&oField7) != OGRERR_NONE || poLayer->CreateField(&oField8) != OGRERR_NONE ||
        poLayer->CreateField(&oField9) != OGRERR_NONE || poLayer->CreateField(&oField10) != OGRERR_NONE || 
        poLayer->CreateField(&oField11) != OGRERR_NONE)
    {
        printf("Creating Name field failed.\n");
        exit(1);
    }
    for (int i = 0; i < tile_infos.size(); ++i) {
        OGRFeature* poFeature;
        poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
        string var_str= to_string(tile_infos[i].ref_var);
        string id_str = to_string(tile_infos[i].tile_id);
        char* cstr = new char[var_str.length() + 1];
        var_str.copy(cstr, var_str.length());
        cstr[var_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        OGR_F_SetFieldString(poFeature, OGR_F_GetFieldIndex(poFeature, "Ref_Var"), cstr);
        delete[] cstr;

        var_str = to_string(tile_infos[i].src_var);
        char* cstr2 = new char[var_str.length() + 1];
        var_str.copy(cstr2, var_str.length());
        cstr2[var_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        OGR_F_SetFieldString(poFeature, OGR_F_GetFieldIndex(poFeature, "Src_Var"), cstr2);
        delete[] cstr2;

        char* cstr1 = new char[id_str.length() + 1];
        id_str.copy(cstr1, id_str.length());
        cstr1[id_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        OGR_F_SetFieldString(poFeature, OGR_F_GetFieldIndex(poFeature, "Tile_ID"), cstr1);
        delete[] cstr1;

        var_str = to_string(tile_infos[i].ref_valid);
        char* cstr3 = new char[var_str.length() + 1];
        var_str.copy(cstr3, var_str.length());
        cstr3[var_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        OGR_F_SetFieldString(poFeature, OGR_F_GetFieldIndex(poFeature, "Ref_Valid"), cstr3);
        delete[] cstr3;

        var_str = to_string(tile_infos[i].src_valid);
        char* cstr4 = new char[var_str.length() + 1];
        var_str.copy(cstr4, var_str.length());
        cstr4[var_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        OGR_F_SetFieldString(poFeature, OGR_F_GetFieldIndex(poFeature, "Src_Valid"), cstr4);
        delete[] cstr4;

        var_str = to_string(tile_infos[i].src_lambda1);
        char* cstr5 = new char[var_str.length() + 1];
        var_str.copy(cstr5, var_str.length());
        cstr5[var_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        OGR_F_SetFieldString(poFeature, OGR_F_GetFieldIndex(poFeature, "Src_Lmd1"), cstr5);
        delete[] cstr5;

        var_str = to_string(tile_infos[i].src_lambda2);
        char* cstr6 = new char[var_str.length() + 1];
        var_str.copy(cstr6, var_str.length());
        cstr6[var_str.length()] = '\0';
        OGR_F_SetFieldString(poFeature, OGR_F_GetFieldIndex(poFeature, "Src_Lmd2"), cstr6);
        delete[] cstr6;

        var_str = to_string(tile_infos[i].src_R);
        char* cstr7 = new char[var_str.length() + 1];
        var_str.copy(cstr7, var_str.length());
        cstr7[var_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        OGR_F_SetFieldString(poFeature, OGR_F_GetFieldIndex(poFeature, "Src_R"), cstr7);
        delete[] cstr7;

        var_str = to_string(tile_infos[i].ref_lambda1);
        char* cstr8 = new char[var_str.length() + 1];
        var_str.copy(cstr8, var_str.length());
        cstr8[var_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        OGR_F_SetFieldString(poFeature, OGR_F_GetFieldIndex(poFeature, "Ref_Lmd1"), cstr8);
        delete[] cstr8;

        var_str = to_string(tile_infos[i].ref_lambda2);
        char* cstr9 = new char[var_str.length() + 1];
        var_str.copy(cstr9, var_str.length());
        cstr9[var_str.length()] = '\0';
        OGR_F_SetFieldString(poFeature, OGR_F_GetFieldIndex(poFeature, "Ref_Lmd2"), cstr9);
        delete[] cstr9;

        var_str = to_string(tile_infos[i].ref_R);
        char* cstr10 = new char[var_str.length() + 1];
        var_str.copy(cstr10, var_str.length());
        cstr10[var_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        OGR_F_SetFieldString(poFeature, OGR_F_GetFieldIndex(poFeature, "Ref_R"), cstr10);
        delete[] cstr10;

        OGRGeometryH geo = OGR_G_CreateGeometry(wkbLinearRing);
        OGR_G_SetPoint_2D(geo, 0, tile_infos[i].utm_bbox[0], tile_infos[i].utm_bbox[2]);
        OGR_G_SetPoint_2D(geo, 1, tile_infos[i].utm_bbox[0], tile_infos[i].utm_bbox[3]);
        OGR_G_SetPoint_2D(geo, 2, tile_infos[i].utm_bbox[1], tile_infos[i].utm_bbox[3]);
        OGR_G_SetPoint_2D(geo, 3, tile_infos[i].utm_bbox[1], tile_infos[i].utm_bbox[2]);
        OGR_G_CloseRings(geo);
        geo = OGR_G_ForceToPolygon(geo);
        OGR_F_SetGeometry(poFeature, geo);
        OGR_G_DestroyGeometry(geo);
        if (poLayer->CreateFeature(poFeature) != OGRERR_NONE)
        {
            printf("Failed to create feature in shapefile.\n");
            exit(1);
        }

        OGRFeature::DestroyFeature(poFeature);
    }

    GDALClose(poDS);

}

void write_tile_info_txt(std::vector<TILE_INFO> tile_infos, string out_path) {
    std::ofstream outfile(out_path);
    for (int  i= 0; i < tile_infos.size(); ++i) {
        outfile << tile_infos[i].tile_id <<"\n";
        
    }
    outfile.close();
}

void write_shift(std::vector<TILE_INFO> tile_infos, string out_path) {
    std::ofstream outfile(out_path);
    for (int i=0;i<tile_infos.size();++i){
        outfile << std::fixed << std::setw(11) << std::setprecision(11) << (double)tile_infos[i].translation[0] << " " <<
            std::fixed << std::setw(11) << std::setprecision(11) << (double)tile_infos[i].translation[1] << " " <<
            std::fixed << std::setw(11) << std::setprecision(11) << (double)tile_infos[i].translation[2] << "\n";
        
    }
    outfile.close();
}

bool cmp_var(pair<float, int>& a, pair<float, int>& b) {
    return a.first > b.first;
}

void ESTIMATE_TRANSFORMATION_CORRS(std::vector<pair<gs::Point, gs::Point>>& corrs, double* rotation, double* translation) {
    gs::Point dynamicMid(0, 0, 0), d_raw, s_raw, d_delta, s_delta;
    gs::Point staticMid(0, 0, 0);
    for (int i = 0; i < corrs.size(); ++i) {
        dynamicMid = dynamicMid + corrs[i].first;
        staticMid = staticMid + corrs[i].second;

    }
    dynamicMid.pos[0] = dynamicMid.pos[0] / corrs.size();
    dynamicMid.pos[1] = dynamicMid.pos[1] / corrs.size();
    dynamicMid.pos[2] = dynamicMid.pos[2] / corrs.size();
    staticMid.pos[0] = staticMid.pos[0] / corrs.size();
    staticMid.pos[1] = staticMid.pos[1] / corrs.size();
    staticMid.pos[2] = staticMid.pos[2] / corrs.size();

    double w[9];
    double U[9]; /// Covariance matrix of Dynamic and Static data.
    double sigma[3];
    double V[9];
    double Diagonal[9];

    double** uSvd = new double* [3];
    double** vSvd = new double* [3];
    uSvd[0] = new double[3];
    uSvd[1] = new double[3];
    uSvd[2] = new double[3];

    vSvd[0] = new double[3];
    vSvd[1] = new double[3];
    vSvd[2] = new double[3];
    gs::clearMatrix(U);
    gs::clearMatrix(V);
    gs::clearMatrix(w);
    for (int i = 0; i < corrs.size(); i++)
    {
        d_raw = corrs[i].first;
        s_raw = corrs[i].second;
        d_delta = d_raw - dynamicMid;
        s_delta = s_raw - staticMid;

        //outerProduct(&d_delta, &s_delta, w);
        gs::outerProduct(&s_delta, &d_delta, w);
        gs::addMatrix(w, U, U);
    }
    gs::copyMatToUV(U, uSvd);
    dsvd(uSvd, 3, 3, sigma, vSvd);
    gs::copyUVtoMat(uSvd, U);
    gs::copyUVtoMat(vSvd, V);

    /// Compute Rotation matrix
    gs::initDiagonal(Diagonal);
    if (gs::determinant_3by3(uSvd) * gs::determinant_3by3(vSvd) < 0) {
        Diagonal[8] = -1;
    }
    //transpose(U);
    gs::transpose(V);
    double tmp[9];
    double tmp_vec3[3];
    gs::clearMatrix(tmp);
    //matrixMult(V, Diagonal,tmp);
    //matrixMult(tmp, U, rotationMatrix);
    gs::matrixMult(U, Diagonal, tmp);
    gs::matrixMult(tmp, V, rotation);

    gs::Point t(0.0, 0.0, 0.0);
    gs::rotate(&dynamicMid, rotation, &t);
    translation[0] = staticMid.pos[0] - t.pos[0];
    translation[1] = staticMid.pos[1] - t.pos[1];
    translation[2] = staticMid.pos[2] - t.pos[2];
    delete[] uSvd[0];
    delete[] uSvd[1];
    delete[] uSvd[2];
    delete[] uSvd;

    delete[] vSvd[0];
    delete[] vSvd[1];
    delete[] vSvd[2];
    delete[] vSvd;
}

void update_utm_bbox(double* utm_bbox, double* rotation, double* translation) {
    utm_bbox[0] += translation[0];
    utm_bbox[1] += translation[0];
    utm_bbox[2] += translation[1];
    utm_bbox[3] += translation[1];
}

void gen_rotation(double* rot, double x_angle,double y_angle,double z_angle) { // in degree
    gs::initDiagonal(rot);
    double rot_tmp1[9];
    double rot_tmp2[9];
    gs::initDiagonal(rot_tmp1);
    gs::initDiagonal(rot_tmp2);

    x_angle = x_angle * 3.1415926 / 180.0;
    y_angle = y_angle * 3.1415926 / 180.0;
    z_angle = z_angle * 3.1415926 / 180.0;
    //Apply x rot
    rot_tmp1[4] = cos(x_angle);
    rot_tmp1[5] = -sin(x_angle);
    rot_tmp1[7] = sin(x_angle);
    rot_tmp1[8] = cos(x_angle);
    gs::matrixMult(rot_tmp1, rot, rot_tmp2);
    //Apply y rot
    gs::initDiagonal(rot_tmp1);
    rot_tmp1[0] = cos(y_angle);
    rot_tmp1[2] = sin(y_angle);
    rot_tmp1[6] = -sin(y_angle);
    rot_tmp1[8] = cos(y_angle);
    gs::matrixMult(rot_tmp1, rot_tmp2, rot);
    //Apply z rot
    gs::initDiagonal(rot_tmp1);
    rot_tmp1[0] = cos(z_angle);
    rot_tmp1[1] = -sin(z_angle);
    rot_tmp1[3] = sin(z_angle);
    rot_tmp1[4] = cos(z_angle);
    gs::matrixMult(rot_tmp1, rot, rot_tmp2);
    std::copy(&rot_tmp2[0], &rot_tmp2[9], &rot[0]);
}
void gen_rotation_zaxis(double* rot, double angle) { // in degree
    gs::initDiagonal(rot);
    angle = angle * 3.1415926 / 180.0;
    rot[0] = cos(angle);
    rot[1] = -sin(angle);
    rot[3] = sin(angle);
    rot[4] = cos(angle);
}

/* [R2,T2]*[R1,T1]=[Rout,Tout]
*/
void merge_two_trans(double* rot1, double* trans1, double* rot2, double* trans2, double* out_rot, double* out_trans) {
	double mid[3];
	gs::matrixMult(rot2, rot1, out_rot);
	gs::rotate_mat(trans1,rot2, mid);
	out_trans[0] = mid[0] + trans2[0];
	out_trans[1] = mid[1] + trans2[1];
	out_trans[2] = mid[2] + trans2[2];
}

void create_dsm(string out_path) {
    int out_width = 70;
    int out_height = 70;
    int min = 45;
    int max = 50;
    GDALAllRegister();
    GDALDriver* pdriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    GDALDataset* out_dataset_ = pdriver->Create(out_path.c_str(), out_width, out_height, 1, GDT_Float64, NULL);
    GDALRasterBand* out_band_ = out_dataset_->GetRasterBand(1);
    double* write_value = new double[out_width*out_height];
    int id = 0;
    for (int row = 0; row < out_height; ++row) {
        for (int col = 0; col < out_width; ++col) {
            id = col + row * out_width;
            if (col < max && col>min && row < max && row>min) {
                write_value[id] = 10.0;
            }
            else {
                write_value[id] = 0.0;
            }
        }
    }
    out_band_->RasterIO(GF_Write, 0, 0, out_width, out_height,
        write_value, out_width, out_height, GDT_Float64, 0, 0);
    GDALClose(out_dataset_);
}

class DSM_REG {
public:
    typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix3X;

    DSM_REG(Parameters par) {
        par_ = par;
        T_rough_ = Eigen::Matrix4d::Identity();
        T_fine_ = Eigen::Matrix4d::Identity();
        T_final_ = Eigen::Matrix4d::Identity();
        gs::initDiagonal(rough_rotation_);
        rough_translation_[0] = 0.0;
        rough_translation_[1] = 0.0;
        rough_translation_[2] = 0.0;
        fine_translation_[0] = 0.0;
        fine_translation_[1] = 0.0;
        fine_translation_[2] = 0.0;

        double src_transform[6]; // [x_origin,x_gsd,_,y_origin,_,-y_gsd]
        double ref_transform[6];
        int xSize, ySize;

        GDALAllRegister();
        pdriver_ = GetGDALDriverManager()->GetDriverByName("GTiff");
        ref_dataset_ = (GDALDataset*)GDALOpen(par_.ref_path_.c_str(), GA_ReadOnly);
        ref_dataset_->GetGeoTransform(ref_transform); 
        xSize = ref_dataset_->GetRasterXSize();
        ySize = ref_dataset_->GetRasterYSize();
        ref_width_ = xSize;
        ref_height_ = ySize;
        ref_utm_bbox_[0] = 0.0;
        ref_utm_bbox_[1] = 0.0 + (xSize + 1) * ref_transform[1];
        ref_utm_bbox_[2] = 0.0 + (ySize - 1) * ref_transform[5];
        ref_utm_bbox_[3] = 0.0;
        ref_reso_ = ref_transform[1];
        GLOBAL_OFFSET_X_m_ = ref_transform[0];
        GLOBAL_OFFSET_Y_m_ = ref_transform[3];


        src_dataset_ = (GDALDataset*)GDALOpen(par_.src_path_.c_str(), GA_ReadOnly);
        src_dataset_->GetGeoTransform(src_transform);
        xSize = src_dataset_->GetRasterXSize();
        ySize = src_dataset_->GetRasterYSize();
        src_width_ = xSize;
        src_height_ = ySize;
        src_utm_bbox_[0] = src_transform[0]- GLOBAL_OFFSET_X_m_;
        src_utm_bbox_[1] = src_transform[0]- GLOBAL_OFFSET_X_m_ + (xSize - 1) * src_transform[1];
        src_utm_bbox_[2] = src_transform[3]- GLOBAL_OFFSET_Y_m_ + (ySize - 1) * src_transform[5];
        src_utm_bbox_[3] = src_transform[3]- GLOBAL_OFFSET_Y_m_;
        src_reso_ = src_transform[1];

        aoi_utm_bbox_[0] = max(ref_utm_bbox_[0], src_utm_bbox_[0]);
        aoi_utm_bbox_[1] = min(ref_utm_bbox_[1], src_utm_bbox_[1]);
        aoi_utm_bbox_[2] = max(ref_utm_bbox_[2], src_utm_bbox_[2]);
        aoi_utm_bbox_[3] = min(ref_utm_bbox_[3], src_utm_bbox_[3]);

        get_aoi_bbox_pixels(aoi_utm_bbox_, ref_utm_bbox_, src_utm_bbox_, ref_width_, ref_height_, src_width_, src_height_, ref_reso_,src_reso_, aoi_ref_bbox_pixel_, aoi_src_bbox_pixel_);

        num_tile_x_ = ceil((aoi_utm_bbox_[1] - aoi_utm_bbox_[0]) / par_.tilesize_m_);
        num_tile_y_ = ceil((aoi_utm_bbox_[3] - aoi_utm_bbox_[2]) / par_.tilesize_m_);

        ref_band_ = ref_dataset_->GetRasterBand(1);
        src_band_ = src_dataset_->GetRasterBand(1);

        par_.src_resolution_x_ = src_transform[1];
        par_.src_resolution_y_ = -src_transform[5];
        par_.ref_resolution_x_ = ref_transform[1];
        par_.ref_resolution_y_ = -ref_transform[5];
        par_.rough_num_pts_ = int(par_.icp_num_pts_ratio_ * src_width_ * src_height_);


        std::cout << "################ Dataset information ################### " << std::endl;
        std::cout << "Global UTM offset (x,y):\n" << GLOBAL_OFFSET_X_m_ << "," << GLOBAL_OFFSET_Y_m_ << std::endl << std::endl;
        std::cout << "src tfw: offset_x:\n" << src_transform[0] << ",offset_y:" << src_transform[3] << ",x_gsd:" << src_transform[1] << ",y_gsd:" << src_transform[5] << std::endl << std::endl;
        std::cout << "src height:" << src_height_ << ",width:" << src_width_ << std::endl << std::endl;
        std::cout << "src offseted bbox (xmin,xmax,ymin,ymax):\n" << src_utm_bbox_[0] << "," << src_utm_bbox_[1] << "," << src_utm_bbox_[2] << "," << src_utm_bbox_[3] << std::endl << std::endl;
        std::cout << "ref tfw: offset_x:\n" << ref_transform[0] << ",offset_y:" << ref_transform[3] << ",x_gsd:" << ref_transform[1] << ",y_gsd:" << ref_transform[5] << std::endl << std::endl;
        std::cout << "ref height:" << ref_height_ << ",width:" << ref_width_ << std::endl << std::endl;
        std::cout << "ref offseted bbox (xmin,xmax,ymin,ymax):\n" << ref_utm_bbox_[0] << "," << ref_utm_bbox_[1] << "," << ref_utm_bbox_[2] << "," << ref_utm_bbox_[3] << std::endl << std::endl;
        std::cout << "aoi offseted bbox (xmin,xmax,ymin,ymax):\n" << aoi_utm_bbox_[0] << "," << aoi_utm_bbox_[1] << "," << aoi_utm_bbox_[2] << "," << aoi_utm_bbox_[3] << std::endl << std::endl;



  //      std::cout << "################ Parameters ################### " << std::endl;
		//std::cout << "ICP point-to-plane: " << par.plane_or_not_ << std::endl;
  //      std::cout << "Rough Align Paras: "<< std::endl;
  //      std::cout << "# Sampling Pts: "<< par_.rough_num_pts_ << std::endl;
  //      std::cout << "NN search half size (pixels): " << par_.search_half_size_pixel_ << std::endl;
  //      std::cout << "Valid Threshold: " << par_.valid_ratio_threshold_rough_ << std::endl;
  //      std::cout << "Rough ICP Max Iters: " << par_.rough_icp_max_iter_ <<'\n'<< std::endl;
  //      std::cout << "Fine Align Paras: " << std::endl;
  //      std::cout << "# Filtered TIles: " << par_.num_filter_tiles_ << std::endl;
  //      std::cout << "Tile size (m): " << par_.tilesize_m_ << std::endl;
  //      std::cout << "Variance Minimal Threshold: " << par_.var_min_ << std::endl;
  //      std::cout << "Valid Threshold: " << par_.valid_ratio_threshold_ <<'\n'<< std::endl;

    };

    Parameters par_;
    GDALDriver* pdriver_;
    GDALDataset* ref_dataset_;
    GDALDataset* src_dataset_;
    GDALRasterBand* ref_band_;
    GDALRasterBand* src_band_;
    int ref_width_, ref_height_, src_width_, src_height_;
    double GLOBAL_OFFSET_X_m_;
    double GLOBAL_OFFSET_Y_m_;
    double ref_utm_bbox_[4];
    double src_utm_bbox_[4];
    double aoi_utm_bbox_[4];
    double ref_reso_, src_reso_;
    int aoi_ref_bbox_pixel_[4];
    int aoi_src_bbox_pixel_[4];
    int num_tile_x_, num_tile_y_;
    std::vector<TILE_INFO> tile_infos_all_, tile_infos_filter_, tile_infos_icp_;

    // for rough align
    Eigen::Matrix4d T_rough_;
    double rough_rotation_[9];
    double rough_translation_[3];
    vector<vector<double>> rough_translations_;


    // for fine align
    Eigen::Matrix4d T_fine_;
    double fine_translation_[3];
    double fine_rmse_;
    double fine_rotation_[9];

    Eigen::Matrix4d T_final_;
    double final_rmse_;
	double final_rot_[9];
	double final_trans_[3];

    // FOR rough align
    std::vector<gs::Point> rough_src_pts_;
    std::vector<std::pair<int, int>> rough_src_index_; // in source image pixel coordinates
    
    static bool cmp_ascend(std::pair<float, int>& a, std::pair<float, int>& b) {
        return a.first < b.first;
    }
    static bool cmp_descend(std::pair<float, int>& a, std::pair<float, int>& b) {
        return a.first > b.first;
    }
    static bool cmp_descend_double(std::pair<double, int>& a, std::pair<double, int>& b) {
        return a.first > b.first;
    }
    static bool cmp_descend_only_double(double& a, double& b) {
        return a > b;
    }
    template <typename T>
    static vector<size_t> sort_indexes(const vector<T>& v) {

        // initialize original index locations
        vector<size_t> idx(v.size());
        iota(idx.begin(), idx.end(), 0);

        // sort indexes based on comparing values in v
        // using std::stable_sort instead of std::sort
        // to avoid unnecessary index re-orderings
        // when v contains elements of equal values 
        stable_sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

        return idx;
    }
    
    void location2pix(double loc_x, double loc_y, double origin_x, double origin_y, double reso_x, double reso_y,int& pix_x,int& pix_y) {
        pix_x = floor((loc_x - origin_x) / reso_x);
        pix_y = floor((loc_y-origin_y) / reso_y);
    }

    void FIND_CORRESPONDENCE(vector<gs::Point> src_pts, Eigen::MatrixXd& corr_src, Eigen::MatrixXd& corr_dst, Eigen::VectorXd& corr_w,
        int search_half_size_pixel, double& RMSE, bool plane_or_not) {
        RMSE = 0;
        bool multi_or_not = false;
        double utm_x, utm_y;
        int grid_x, grid_y;
        int search_xmin, search_xmax, search_ymax, search_ymin;
        int search_width, search_height;
        double* ref_search_area = new double[(2 * search_half_size_pixel + 1) * (2 * search_half_size_pixel + 1)];
        double src_pts_num = src_pts.size();
        std::vector<CORR> corrs;
        for (int i = 0; i < src_pts.size(); ++i) {
            double src_pt[3] = { src_pts[i].pos[0],src_pts[i].pos[1],src_pts[i].pos[2] };
            //find pixel location in ref_image coordinate
            location2pix(src_pt[0], src_pt[1], ref_utm_bbox_[0], ref_utm_bbox_[3], par_.ref_resolution_x_, -par_.ref_resolution_y_, grid_x, grid_y);
            if (grid_x < 0 || grid_x >= ref_width_ || grid_y < 0 || grid_y >= ref_height_) {
                continue;
            }
            search_xmin = grid_x - search_half_size_pixel;
            search_xmax = grid_x + search_half_size_pixel;
            search_ymin = grid_y - search_half_size_pixel;
            search_ymax = grid_y + search_half_size_pixel;
            search_xmin = (search_xmin < 0) ? 0 : search_xmin;
            search_xmax = (search_xmax >= ref_width_) ? ref_width_ - 1 : search_xmax;
            search_ymin = (search_ymin < 0) ? 0 : search_ymin;
            search_ymax = (search_ymax >= ref_height_) ? ref_height_ - 1 : search_ymax;
            search_width = search_xmax - search_xmin + 1;
            search_height = search_ymax - search_ymin + 1;

            ref_band_->RasterIO(GF_Read, search_xmin, search_ymin,
                search_width, search_height,
                ref_search_area, search_width, search_height, GDT_Float64, 0, 0);

            //construct kd-tree
            PointCloud<double> cloud;
            for (int row = 0; row < search_height ; ++row) {
                for (int col = 0; col < search_width; ++col) {
                    int idx = col + row * search_width;
                    if (isnan(ref_search_area[idx])) {
                        continue;
                    }
                    double x = ref_utm_bbox_[0] + (search_xmin + col+0.5) * par_.ref_resolution_x_;
                    double y = ref_utm_bbox_[3] - (search_ymin + row+0.5) * par_.ref_resolution_y_;
                    double z = ref_search_area[idx];
                    PointCloud<double>::Point p;
                    p.x = x, p.y = y, p.z = z;
                    cloud.pts.push_back(p);
                }
            }
            if (cloud.pts.size() == 0) {
                continue;
            }
            using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
                nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
                PointCloud<double>, 3 /* dim */>;
            my_kd_tree_t kd_index(3 /*dim*/, cloud, { 10 /* max leaf */ });
            const size_t                   num_results = 1;
            size_t                         ret_index;
            double                          out_dist_sqr;
            nanoflann::KNNResultSet<double> resultSet(num_results);
            resultSet.init(&ret_index, &out_dist_sqr);
            kd_index.findNeighbors(resultSet, &src_pt[0],nanoflann::SearchParams(10));

            gs::Point ref_pt(cloud.pts[ret_index].x, cloud.pts[ret_index].y, cloud.pts[ret_index].z);
            gs::Point src_pt_point(src_pt[0], src_pt[1], src_pt[2]);
            CORR corr;
            corr.ref = ref_pt;
            corr.src = src_pt_point;
            corr.w = 1;
            corrs.push_back(corr);
        }
        //reweight the corr
        for (int i = 0; i < corrs.size(); ++i) {
            corrs[i].w /= double(corrs.size());
        }
        // reject outlier
        if (!multi_or_not) {
            CORR_REJECTOR1(corrs);
        }

        //compute rmse
        if (plane_or_not) {
            for (int i = 0; i < corrs.size(); ++i) {
                RMSE += corrs[i].w * (
                    pow((corrs[i].src.pos[0] - corrs[i].ref.pos[0]) * corrs[i].ref_norm[0], 2) +
                    pow((corrs[i].src.pos[1] - corrs[i].ref.pos[1]) * corrs[i].ref_norm[1], 2) +
                    pow((corrs[i].src.pos[2] - corrs[i].ref.pos[2]) * corrs[i].ref_norm[2], 2)
                    );
            }
            RMSE = sqrt(RMSE);
        }
        else {
            for (int i = 0; i < corrs.size(); ++i) {
                CORR corr = corrs[i];
                RMSE += corrs[i].w * (pow(corrs[i].src.pos[0] - corrs[i].ref.pos[0], 2.0) +
                    pow(corrs[i].src.pos[1] - corrs[i].ref.pos[1], 2.0) +
                    pow(corrs[i].src.pos[2] - corrs[i].ref.pos[2], 2.0));
                if (isnan(RMSE)) {
                    int a = 1;
                }
            }
            RMSE = sqrt(RMSE);
        }
        if (corrs.size() == 0) {
            RMSE = 9999;
        }

        //transform to eigen type
        corr_src.resize(3, corrs.size());
        corr_dst.resize(3, corrs.size());
        corr_w.resize(corrs.size());
        for (int i = 0; i < corrs.size(); ++i) {
            corr_src(0, i) = corrs[i].src.pos[0];
            corr_src(1, i) = corrs[i].src.pos[1];
            corr_src(2, i) = corrs[i].src.pos[2];
            corr_dst(0, i) = corrs[i].ref.pos[0];
            corr_dst(1, i) = corrs[i].ref.pos[1];
            corr_dst(2, i) = corrs[i].ref.pos[2];
            corr_w(i) = corrs[i].w;
        }
        delete[] ref_search_area;
    }

    void FIND_CORRESPONDENCE_ADAPTIVE(vector<gs::Point> src_pts, Eigen::MatrixXd& corr_src, 
                                        Eigen::MatrixXd& corr_dst, Eigen::VectorXd& corr_w,
                                        double& RMSE, bool plane_or_not) {
        RMSE = 0;
        bool robust_reject = false;
        double utm_x, utm_y;
        int grid_x, grid_y;
        int search_xmin, search_xmax, search_ymax, search_ymin;
        int search_width, search_height;
        double* ref_search_pt = new double[1];
        double src_pts_num = src_pts.size();
        std::vector<CORR> corrs;
        if (par_.outlier_thresh_ != -1) {
            robust_reject = true;
        }
        for (int i = 0; i < src_pts.size(); ++i) {
            double src_pt[3] = { src_pts[i].pos[0],src_pts[i].pos[1],src_pts[i].pos[2] };
            //find pixel location in ref_image coordinate
            location2pix(src_pt[0], src_pt[1], ref_utm_bbox_[0], ref_utm_bbox_[3], par_.ref_resolution_x_, -par_.ref_resolution_y_, grid_x, grid_y);
            if (grid_x < 0 || grid_x >= ref_width_ || grid_y < 0 || grid_y >= ref_height_) {
                continue;
            }

            ref_band_->RasterIO(GF_Read, grid_x, grid_y,1, 1, ref_search_pt, 1, 1, GDT_Float64, 0, 0);
            if (isnan(ref_search_pt[0]) || ref_search_pt[0]==-9999) {
                continue;
            }
            int search_half_size_pixel_x = ceil(abs(src_pt[2] - ref_search_pt[0])/par_.ref_resolution_x_);
            int search_half_size_pixel_y = ceil(abs(src_pt[2] - ref_search_pt[0])/par_.ref_resolution_y_);
            if (search_half_size_pixel_x == 1 && search_half_size_pixel_y == 1) {
                double x = ref_utm_bbox_[0] + (grid_x + 0.5) * par_.ref_resolution_x_;
                double y = ref_utm_bbox_[3] - (grid_y + 0.5) * par_.ref_resolution_y_;
                double z = ref_search_pt[0];
                gs::Point ref_pt(x,y,z);
                gs::Point src_pt_point(src_pt[0], src_pt[1], src_pt[2]);
                CORR corr;
                corr.ref = ref_pt;
                corr.src = src_pt_point;
                corr.w = 1;
                corrs.push_back(corr);
            }
            else
            {
                double* ref_search_area = new double[(2 * search_half_size_pixel_x + 1) * (2 * search_half_size_pixel_y + 1)];
                search_xmin = grid_x - search_half_size_pixel_x;
                search_xmax = grid_x + search_half_size_pixel_x;
                search_ymin = grid_y - search_half_size_pixel_y;
                search_ymax = grid_y + search_half_size_pixel_y;
                search_xmin = (search_xmin < 0) ? 0 : search_xmin;
                search_xmax = (search_xmax >= ref_width_) ? ref_width_ - 1 : search_xmax;
                search_ymin = (search_ymin < 0) ? 0 : search_ymin;
                search_ymax = (search_ymax >= ref_height_) ? ref_height_ - 1 : search_ymax;
                search_width = search_xmax - search_xmin + 1;
                search_height = search_ymax - search_ymin + 1;

                ref_band_->RasterIO(GF_Read, search_xmin, search_ymin,
                    search_width, search_height,
                    ref_search_area, search_width, search_height, GDT_Float64, 0, 0);

                //construct kd-tree
                //PointCloud<double> cloud;
                double min_dis = 0;
                int min_idx = 0;
                double min_x = 0 ;
                double min_y = 0;
                double min_z = 0 ;
                for (int row = 0; row < search_height; ++row) {
                    for (int col = 0; col < search_width; ++col) {
                        int idx = col + row * search_width;
                        if (isnan(ref_search_area[idx]) || ref_search_area[idx]==-9999 ) {
                            continue;
                        }
                        double x = ref_utm_bbox_[0] + (search_xmin + col + 0.5) * par_.ref_resolution_x_;
                        double y = ref_utm_bbox_[3] - (search_ymin + row + 0.5) * par_.ref_resolution_y_;
                        double z = ref_search_area[idx];
                        double dis = (x - src_pt[0]) * (x - src_pt[0]) + (y - src_pt[1]) * (y - src_pt[1]) + (z - src_pt[2]) * (z - src_pt[2]);
                        if (idx == 0) {
                            min_dis = dis;
                            min_x = x;
                            min_y = y;
                            min_z = z;
                        }
                        else {
                            if (dis < min_dis) {
                                min_dis = dis;
                                min_idx = idx;
                                min_x = x;
                                min_y = y;
                                min_z = z;
                            }
                        }
                        //PointCloud<double>::Point p;
                        //p.x = x, p.y = y, p.z = z;
                        //cloud.pts.push_back(p);
                    }
                }
                //if (cloud.pts.size() == 0) {
                //    continue;
                //}
                //using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
                //    nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
                //    PointCloud<double>, 3 /* dim */>;
                //my_kd_tree_t kd_index(3 /*dim*/, cloud, { 10 /* max leaf */ });
                //const size_t                   num_results = 1;
                //size_t                         ret_index;
                //double                          out_dist_sqr;
                //nanoflann::KNNResultSet<double> resultSet(num_results);
                //resultSet.init(&ret_index, &out_dist_sqr);
                //kd_index.findNeighbors(resultSet, &src_pt[0], nanoflann::SearchParams(10));

                //gs::Point ref_pt(cloud.pts[ret_index].x, cloud.pts[ret_index].y, cloud.pts[ret_index].z);
                gs::Point ref_pt(min_x, min_y, min_z);
                gs::Point src_pt_point(src_pt[0], src_pt[1], src_pt[2]);
                CORR corr;
                corr.ref = ref_pt;
                corr.src = src_pt_point;
                corr.w = 1;
                corrs.push_back(corr);
                delete[] ref_search_area;
            }
        }
        // robust reject outlier correspondences
        if (par_.outlier_thresh_ != -1) {
            CORR_REJECTOR1(corrs);
        }
        //reweight the corr
        for (int i = 0; i < corrs.size(); ++i) {
            corrs[i].w /= double(corrs.size());
        }

        //compute rmse
        if (plane_or_not) {
            for (int i = 0; i < corrs.size(); ++i) {
                RMSE += corrs[i].w * (
                    pow((corrs[i].src.pos[0] - corrs[i].ref.pos[0]) * corrs[i].ref_norm[0], 2) +
                    pow((corrs[i].src.pos[1] - corrs[i].ref.pos[1]) * corrs[i].ref_norm[1], 2) +
                    pow((corrs[i].src.pos[2] - corrs[i].ref.pos[2]) * corrs[i].ref_norm[2], 2)
                    );
            }
            RMSE = sqrt(RMSE);
        }
        else {
            for (int i = 0; i < corrs.size(); ++i) {
                CORR corr = corrs[i];
                RMSE += (pow(corrs[i].src.pos[0] - corrs[i].ref.pos[0], 2.0) +
                    pow(corrs[i].src.pos[1] - corrs[i].ref.pos[1], 2.0) +
                    pow(corrs[i].src.pos[2] - corrs[i].ref.pos[2], 2.0));

            }
            RMSE /= corrs.size();
            RMSE = sqrt(RMSE);
        }
        if (corrs.size() == 0) {
            RMSE = -1;
        }

        //transform to eigen type
        corr_src.resize(3, corrs.size());
        corr_dst.resize(3, corrs.size());
        corr_w.resize(corrs.size());
        for (int i = 0; i < corrs.size(); ++i) {
            corr_src(0, i) = corrs[i].src.pos[0];
            corr_src(1, i) = corrs[i].src.pos[1];
            corr_src(2, i) = corrs[i].src.pos[2];
            corr_dst(0, i) = corrs[i].ref.pos[0];
            corr_dst(1, i) = corrs[i].ref.pos[1];
            corr_dst(2, i) = corrs[i].ref.pos[2];
            corr_w(i) = corrs[i].w;
        }
        delete[] ref_search_pt;

    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="X">source points, 3xN</param>
    /// <param name="Y">target points, 3xN</param>
    /// <param name="W">weight, N</param>
    /// <param name="R">estimate rotation, 3x3</param>
    /// <param name="T">estimate translation, 3</param>
    void ESTIMATE_TRANSFORMATION_WEIGHTED_CORRS(Eigen::MatrixXd& X, Eigen::MatrixXd& Y, Eigen::VectorXd& W, 
                                                Eigen::Matrix4d& trans,std::string type) {
        int dim = X.rows();
        /// Normalize weight vector
        Eigen::VectorXd w_normalized = W / W.sum();
        /// De-mean
        Eigen::VectorXd X_mean(dim), Y_mean(dim);
        for (int i = 0; i < dim; ++i) {
            X_mean(i) = (X.row(i).array() * w_normalized.transpose().array()).sum();
            Y_mean(i) = (Y.row(i).array() * w_normalized.transpose().array()).sum();
        }

        if (type == "translation") {
            trans.block(0, 3, 3, 1) = Y_mean - X_mean;
        }
        else if(type == "rigid") {
            X.colwise() -= X_mean;
            Y.colwise() -= Y_mean;
            /// Compute transformation
            Eigen::MatrixXd sigma = X * w_normalized.asDiagonal() * Y.transpose();
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(sigma, Eigen::ComputeFullU | Eigen::ComputeFullV);
            if (svd.matrixU().determinant() * svd.matrixV().determinant() < 0.0) {
                Eigen::VectorXd S = Eigen::VectorXd::Ones(dim); S(dim - 1) = -1.0;
                trans.block(0,0,3,3) = svd.matrixV() * S.asDiagonal() * svd.matrixU().transpose();
            }
            else {
                trans.block(0, 0, 3, 3) = svd.matrixV() * svd.matrixU().transpose();
            }
            trans.block(0, 3, 3, 1) = Y_mean - trans.block(0, 0, 3, 3) * X_mean;
            /// Re-apply mean
            X.colwise() += X_mean;
            Y.colwise() += Y_mean;
        }


    }

    void CORR_REJECTOR(std::vector<pair<gs::Point, gs::Point>>& corrs) {
        int num_trimmed = (int)corrs.size() * par_.trimmed_ratio_;
        std::vector<pair<gs::Point, gs::Point>> corrs_tmp;
        std::vector<std::pair<float, int>> corr_dis_index;
        float dis;
        for (int i = 0; i < corrs.size(); ++i) {
            dis = std::sqrt(pow(corrs[i].first.pos[0] - corrs[i].second.pos[0], 2.0) +
                pow(corrs[i].first.pos[1] - corrs[i].second.pos[1], 2.0) +
                pow(corrs[i].first.pos[2] - corrs[i].second.pos[2], 2.0));
            corr_dis_index.push_back(std::make_pair(dis, i));
        }
        std::sort(corr_dis_index.begin(), corr_dis_index.end(), cmp_ascend);
        for (int i = 0; i < num_trimmed; ++i) {
            corrs_tmp.push_back(corrs[corr_dis_index[i].second]);
        }
        corrs.clear();
        corrs = corrs_tmp;
    }
    void CORR_REJECTOR1(std::vector<CORR>& corrs) {
        std::vector<CORR> corrs_tmp;
        std::vector<std::pair<float, int>> corr_dis_index;
        float dis;
        for (int i = 0; i < corrs.size(); ++i) {
            dis = std::sqrt(pow(corrs[i].src.pos[0] - corrs[i].ref.pos[0], 2.0) +
                pow(corrs[i].src.pos[1] - corrs[i].ref.pos[1], 2.0) +
                pow(corrs[i].src.pos[2] - corrs[i].ref.pos[2], 2.0));
            if (dis > par_.outlier_thresh_) {
                continue;
            }
            corrs_tmp.push_back(corrs[i]);
        }
        corrs.clear();
        corrs = corrs_tmp;
    }
    
    void START(std::string debug_dir) {
        print_info();
        ADAPTIVE_ALIGN(debug_dir);
        STEP_END();
    }

    void ADAPTIVE_ALIGN(string dir_path) {
        // Sample src pts, uniformly distributed in the aoi region
        double margin_ratio = 0.1;
        int aoi_src_bbox_pixel_margin[4];
        int margin_width = (aoi_src_bbox_pixel_[1] - aoi_src_bbox_pixel_[0]) * margin_ratio / 2;
        int margin_height = (aoi_src_bbox_pixel_[2] - aoi_src_bbox_pixel_[3]) * margin_ratio / 2;
        aoi_src_bbox_pixel_margin[0] = aoi_src_bbox_pixel_[0] + margin_width;
        aoi_src_bbox_pixel_margin[1] = aoi_src_bbox_pixel_[1] - margin_width;
        aoi_src_bbox_pixel_margin[2] = aoi_src_bbox_pixel_[2] - margin_height;
        aoi_src_bbox_pixel_margin[3] = aoi_src_bbox_pixel_[3] + margin_height;

        double widht_height_ratio = (aoi_utm_bbox_[1] - aoi_utm_bbox_[0]) / (aoi_utm_bbox_[3] - aoi_utm_bbox_[2]);
        int grid_height = ceil(std::sqrt((double)par_.rough_num_pts_ / widht_height_ratio));
        int grid_width = (int)((double)grid_height * widht_height_ratio);
        int grid_step_x_pixel = (aoi_src_bbox_pixel_margin[1] - aoi_src_bbox_pixel_margin[0]) / grid_width;
        int grid_step_y_pixel = (aoi_src_bbox_pixel_margin[2] - aoi_src_bbox_pixel_margin[3]) / grid_height;
        int grid_x, grid_y;
        double utm_x, utm_y;
        double* pts_z = new double[1];
        for (int row = 0; row < grid_height; ++row) {
            for (int col = 0; col < grid_width; ++col) {
                grid_x = aoi_src_bbox_pixel_margin[0] + col * grid_step_x_pixel;
                grid_y = aoi_src_bbox_pixel_margin[3] + row * grid_step_y_pixel;
                utm_x = src_utm_bbox_[0] + (grid_x+0.5) * par_.src_resolution_x_;
                utm_y = src_utm_bbox_[3] - (grid_y+0.5) * par_.src_resolution_y_;
                src_band_->RasterIO(GF_Read, grid_x, grid_y, 1, 1, pts_z, 1, 1, GDT_Float64, 0, 0);
                if (isnan(*pts_z)) {
                    continue;
                }
                rough_src_index_.push_back(make_pair(grid_x, grid_y));
                gs::Point src_pt(utm_x, utm_y, *pts_z);
                rough_src_pts_.push_back(src_pt);
            }
        }
        delete[] pts_z;
        std::cout << "#pts :"<<rough_src_pts_.size()<<"used in ICP";
        // ICP
        string rough_icp_log = dir_path + "rough_icp.log";
        string rough_trans_txt = dir_path + "rough_icp_trans.txt";
        double rough_rmse;
        gs::initDiagonal(rough_rotation_);
        rough_translation_[0] = 0;
        rough_translation_[1] = 0;
        rough_translation_[2] = 0;
        T_final_ = Eigen::Matrix4d::Identity();
        ADAPTIVE_ICP(rough_src_pts_, T_final_, par_.type_, par_.search_half_size_pixel_, par_.rough_icp_max_iter_, par_.rough_icp_rmse_threshold_,
            par_.plane_or_not_, final_rmse_, rough_icp_log);
    }

    void STEP_END() {
        GDALClose(ref_dataset_);
        GDALClose(src_dataset_);
    }
    
    // correspondence is based on same x,y coordinates
    void ADAPTIVE_ICP(std::vector<gs::Point> src_pts, Eigen::Matrix4d& T, std::string type, int init_search_half_pixels, int max_iter, double rmse_threshold, double plane_or_not, double& final_rmse, string debug_path) {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        update_ptcloud(src_pts, T, 0, 0, 0, 0);
        Eigen::MatrixXd X, Y;
        Eigen::VectorXd W;
        double corr_rmse = 0;
        FIND_CORRESPONDENCE_ADAPTIVE(src_pts, X, Y, W, corr_rmse, plane_or_not);
        final_rmse = corr_rmse;
        //Estmate Transformation
        ofstream ofs;
        if (debug_path != "" && par_.verbose_ > 0) {
            ofs.open(debug_path);
        }
        for (int iter = 0; iter < max_iter; ++iter) {
            std::cout << std::setprecision(11) << "ICP #Iteration: " << iter << "# corrs: " << X.cols() << " RMSE: " << corr_rmse << "Trans: " << T(0, 3) << "," << T(1, 3) << "," << T(2, 3) << std::endl;
            if (debug_path != "" && par_.verbose_ > 0) {
                ofs << std::setprecision(11) << corr_rmse << " " << T(0, 3) << " " << T(1, 3) << " " << T(2, 3) << std::endl;
            }
            Eigen::Matrix4d transform;
            double rotation[9];
            double translation[3];
            ESTIMATE_TRANSFORMATION_WEIGHTED_CORRS(X, Y, W, transform, type);

            T.block(0, 0, 3, 3) = transform.block(0, 0, 3, 3) * T.block(0, 0, 3, 3);
            T.block(0, 3, 3, 1) = transform.block(0, 0, 3, 3) * T.block(0, 3, 3, 1) + transform.block(0, 3, 3, 1);

            // Update source pts
            update_ptcloud(src_pts, transform, 0.0, 0.0, 0.0, 0.0);

            // FInd correspondence
            double pre_rmse = corr_rmse;
            FIND_CORRESPONDENCE_ADAPTIVE(src_pts, X, Y, W, corr_rmse, plane_or_not);

            if (abs(pre_rmse - corr_rmse) < rmse_threshold) {
                final_rmse = (corr_rmse < final_rmse) ? corr_rmse : final_rmse;
                break;
            }
            final_rmse = (corr_rmse < final_rmse) ? corr_rmse : final_rmse;
        }
        if (debug_path != "" && par_.verbose_ > 0) {
            ofs.close();
        }

        //print
        std::cout << "Final RMSE (< " << par_.outlier_thresh_ << "m): " << final_rmse << std::endl;
        std::cout << "ICP result rotation: \n" << T.block(0, 0, 3, 3) << "\n" <<
            "ICP result translation\n" << T.block(0, 3, 3, 1) << endl;
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds> (end - begin).count()/1000000.0 << "[s]" << std::endl;

    }

    void print_info() {
        std::cout << "Registration Start\n";
        if (par_.multi_or_not_) {
            std::cout << "Multi-link function on\n";
        }
        else {
            std::cout << "Multi-link function off\n";
        }
        if (par_.plane_or_not_) {
            std::cout << "Point-to-Plane function on\n";
        }
        else {
            std::cout << "Point-to-Plane function off\n";
        }
        if (par_.brute_force_search_or_not_) {
            std::cout << "Brute-Force search on X-Y plane on\n";
        }
        else {
            std::cout << "Brute-Force search on X-Y plane off\n";
        }

    }
};

class DSM_TRANSFORM {
public:
    DSM_TRANSFORM(std::string in_path,std::string out_path,Eigen::Matrix4d Transform,double offx,double offy) {
        transform_ = Transform;
        in_zmin_ = -9999;
        in_zmax_ = -9999;
        int tile_size = 3000;
        //std::copy(&rotation[0], &rotation[9], &rotation_[0]);
        //std::copy(&translation[0], &translation[3], &translation_[0]);
        GLOBAL_OFFSET_X_m_ = offx;
        GLOBAL_OFFSET_Y_m_ = offy;
        double transform[6],transform_out[6];
        int xSize, ySize;
        out_path_ = out_path;
        in_path_ = in_path;
        GDALAllRegister();
        pdriver_ = GetGDALDriverManager()->GetDriverByName("GTiff");
        in_dataset_ = (GDALDataset*)GDALOpen(in_path_.c_str(), GA_ReadOnly);
        in_dataset_->GetGeoTransform(transform);
        xSize = in_dataset_->GetRasterXSize();
        ySize = in_dataset_->GetRasterYSize();
        in_width_ = xSize;
        in_height_ = ySize;
        in_utm_bbox_[0] = transform[0];
        in_utm_bbox_[1] = transform[0] + xSize * transform[1];
        in_utm_bbox_[2] = transform[3] + ySize * transform[5];
        in_utm_bbox_[3] = transform[3];
        reso_x_ = transform[1];
        reso_y_ = -transform[5];
        in_datatype_=GDALGetRasterDataType(GDALGetRasterBand(in_dataset_, 1));
        in_band_ = in_dataset_->GetRasterBand(1);
        //Collect zmin and zmax of in_raster
        int num_tile_x = ceil((double)in_width_ / (double)tile_size);
        int num_tile_y = ceil((double)in_height_ / (double)tile_size);
        int fill_x = 0;
        int fill_y = 0;
        int init_cnt = 0;
        for (int row = 0; row < num_tile_y; ++row) {
            for (int col = 0; col < num_tile_x; ++col) {
                
                fill_x = (tile_size * (col + 1) <= in_width_) ? tile_size * (col + 1) : in_width_;
                fill_x -= tile_size * col;
                fill_y = (tile_size * (row + 1) <= in_height_) ? tile_size * (row + 1) : in_height_;
                fill_y -= tile_size * row;
                double* in_data = new double[fill_x*fill_y];
                in_band_->RasterIO(GF_Read, col, row, fill_x, fill_y,
                    in_data, fill_x, fill_y, GDT_Float64, 0, 0);
                for (int i = 0; i < fill_x * fill_y; ++i) {
                    if (isnan(in_data[i])) {
                        continue;
                    }
                    init_cnt++;
                    if (init_cnt==1) {
                        in_zmin_ = in_data[i];
                        in_zmax_ = in_data[i];
                    }
                    if (in_data[i] > in_zmax_) {
                        in_zmax_ = in_data[i];
                    }
                    if (in_data[i] < in_zmin_) {
                        in_zmin_ = in_data[i];
                    }
                }
                delete[] in_data;
            }
        }
        // compute output bbox size
        //double pt1[3], pt2[3], pt3[3], pt4[3], pt5[3], pt6[3], pt7[3], pt8[3],
        //       pt1_trans[3],pt2_trans[3],pt3_trans[3],pt4_trans[3], pt5_trans[3], pt6_trans[3], pt7_trans[3], pt8_trans[3];
        Eigen::Vector3d pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8,
            pt1_trans, pt2_trans, pt3_trans, pt4_trans, pt5_trans, pt6_trans, pt7_trans, pt8_trans;
        pt1 << in_utm_bbox_[0],  in_utm_bbox_[3],in_zmin_;
        pt5<< in_utm_bbox_[0],in_utm_bbox_[3],in_zmax_;
        pt2<<in_utm_bbox_[1],in_utm_bbox_[3], in_zmin_;
        pt6<<in_utm_bbox_[1],in_utm_bbox_[3], in_zmax_;
        pt3<<in_utm_bbox_[1],in_utm_bbox_[2],in_zmin_;
        pt7<<in_utm_bbox_[1],in_utm_bbox_[2],in_zmax_;
        pt4<<in_utm_bbox_[0], in_utm_bbox_[2],in_zmin_;
        pt8<<in_utm_bbox_[0],in_utm_bbox_[2],in_zmax_;
        point_trans(pt1, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, pt1_trans);
        point_trans(pt2, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, pt2_trans);
        point_trans(pt3, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, pt3_trans);
        point_trans(pt4, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, pt4_trans);
        point_trans(pt5, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, pt5_trans);
        point_trans(pt6, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, pt6_trans);
        point_trans(pt7, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, pt7_trans);
        point_trans(pt8, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, pt8_trans);
        out_utm_bbox_[0] = min4(pt1_trans(0), pt4_trans(0), pt5_trans(0), pt8_trans(0));
        out_utm_bbox_[1] = max4(pt2_trans(0), pt3_trans(0), pt6_trans(0), pt7_trans(0));
        out_utm_bbox_[2] = min4(pt3_trans(1), pt4_trans(1), pt7_trans(1), pt8_trans(1));
        out_utm_bbox_[3] = max4(pt1_trans(1), pt2_trans(1), pt5_trans(1), pt6_trans(1));
        out_width_ = floor((out_utm_bbox_[1] - out_utm_bbox_[0]+transform[1]/2) / transform[1])+1;
        out_height_ = floor((out_utm_bbox_[3] - out_utm_bbox_[2] + transform[1] / 2) / transform[1]) + 1;

        std::copy(std::begin(transform), std::end(transform), std::begin(transform_out));
        transform_out[0] = out_utm_bbox_[0];
        transform_out[3] = out_utm_bbox_[3];
        // Create Empty output dsm
        out_dataset_ = pdriver_->Create(out_path_.c_str(), out_width_, out_height_, 1, in_datatype_, NULL);
        out_dataset_->SetGeoTransform(transform_out);
        out_band_ = out_dataset_->GetRasterBand(1);
        num_tile_x = ceil((double)out_width_ / (double)tile_size);
        num_tile_y= ceil((double)out_height_ / (double)tile_size);
        for (int row = 0; row < num_tile_y; ++row) {
            for (int col = 0; col < num_tile_x; ++col) {
                fill_x = (tile_size * (col + 1) <=out_width_)? tile_size * (col + 1): out_width_;
                fill_x -= tile_size * col;
                fill_y = (tile_size * (row + 1) <= out_height_) ? tile_size * (row + 1) : out_height_;
                fill_y -= tile_size * row;
                rasterio_write_nan(out_band_, tile_size * col, tile_size * row, fill_x, fill_y);
            }
        }
        std::string out_tfw_path = out_path_;
        strutil::replace_last(out_tfw_path, ".tif", ".tfw");
        export_tfw(out_tfw_path, transform_out);
        //GDALClose(out_dataset_);
    };
    std::string in_path_, out_path_;
    GDALDriver* pdriver_;
    GDALDataset* out_dataset_;
    GDALDataset* in_dataset_;
    GDALRasterBand* in_band_;
    GDALRasterBand* out_band_;
    GDALDataType in_datatype_;
    int out_width_, out_height_, in_width_, in_height_;
    double GLOBAL_OFFSET_X_m_;
    double GLOBAL_OFFSET_Y_m_;
    double reso_x_;
    double reso_y_;
    double out_utm_bbox_[4];
    double in_utm_bbox_[4];
    double in_zmin_, in_zmax_;
    Eigen::Matrix4d transform_;

    void START() {
        double* write_value = new double[1];
        double* read_value = new double[1];
        Eigen::Vector3d pt, pt_out;
        int col_out, row_out;
        for (int row = 0; row < in_height_; ++row)
        {
            for (int col = 0; col < in_width_; ++col) {
                pt[0] = in_utm_bbox_[0] + reso_x_ * (col+0.5);
                pt[1] = in_utm_bbox_[3] - reso_y_ * (row+0.5);
                in_band_->RasterIO(GF_Read, col, row, 1, 1,
                    read_value, 1, 1, GDT_Float64, 0, 0);
                pt[2] = read_value[0];
                if (isnan(pt[2])) {
                    continue;
                }
                point_trans(pt, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, pt_out);
                col_out = floor((pt_out[0] - out_utm_bbox_[0]) / reso_x_);
                row_out= floor((out_utm_bbox_[3]- pt_out[1]) / reso_y_);
                if (row_out < 0) {
                    int a = 0;
                }
                //std::cout << col_out << "," << row_out << std::endl;
                out_band_->RasterIO(GF_Read, col_out, row_out, 1, 1,
                    read_value, 1, 1, GDT_Float64, 0, 0);
                if (pt_out[2] > read_value[0] || isnan(read_value[0])) {
                    out_band_->RasterIO(GF_Write, col_out, row_out, 1, 1,
                        &pt_out[2], 1, 1, GDT_Float64, 0, 0);
                }
            }
        }
        GDALClose(in_dataset_);
        GDALClose(out_dataset_);
    }
    double min4(double a, double b, double c, double d) {
        double tmp;
        tmp = min(a, b);
        tmp = min(tmp, c);
        tmp = min(tmp, d);
        return tmp;
    }
    double max4(double a, double b, double c, double d) {
        double tmp;
        tmp = max(a, b);
        tmp = max(tmp, c);
        tmp = max(tmp, d);
        return tmp;
    }
    void gen_rotation_zaxis(double* rot, double angle) { // in degree
        gs::initDiagonal(rot);
        angle = angle * 3.1415926 / 180.0;
        rot[0] = cos(angle);
        rot[1] = -sin(angle);
        rot[3] = sin(angle);
        rot[4] = cos(angle);
    }
    void rasterio_write_nan(GDALRasterBand* band, int offsetx, int offsety, int xsize, int ysize) {
        int cnt = xsize * ysize;
        GDALDataType dt_in = GDALGetRasterDataType(band);
        int bits = GDALGetDataTypeSize(dt_in);
        //cout << "rasterio_write_nan: " << bits << endl;
        switch (bits) {
        case 16: {
            //GInt16* buf16 = new GInt16[cnt];
            GInt16* buf16 = (GInt16*)CPLMalloc(sizeof(GInt16) * cnt);
            for (int i = 0; i < cnt; ++i) {
                buf16[i] = std::numeric_limits<GInt16>::quiet_NaN();
            }
            band->RasterIO(GF_Write, offsetx, offsety, xsize, ysize, buf16, xsize, ysize, GDT_Int16, 0, 0);
            //delete[] buf16;
            CPLFree(buf16);
            break;
        }
        case 32: {
            //float* buf32 = new float[cnt];
            float* buf32 = (float*)CPLMalloc(sizeof(float) * cnt);
            for (int i = 0; i < cnt; ++i) {
                buf32[i] = std::numeric_limits<float>::quiet_NaN();
            }
            band->RasterIO(GF_Write, offsetx, offsety, xsize, ysize, buf32, xsize, ysize, GDT_Float32, 0, 0);
            CPLFree(buf32);
            //delete[] buf32;
            break;

        }

        case 64: {
            float* buf64 = new float[cnt];
            for (int i = 0; i < cnt; ++i) {
                buf64[i] = std::numeric_limits<double>::quiet_NaN();
            }
            band->RasterIO(GF_Write, offsetx, offsety, xsize, ysize, buf64, xsize, ysize, GDT_Float64, 0, 0);
            break;
        }

        }

    }
    void export_tfw(string path, double* transform) {
        std::ofstream out(path.c_str());
        out << std::fixed << std::setprecision(8) << transform[1] << endl;
        out << std::fixed << std::setprecision(8) << transform[2] << endl;
        out << std::fixed << std::setprecision(8) << transform[4] << endl;
        out << std::fixed << std::setprecision(8) << transform[5] << endl;
        out << std::fixed << std::setprecision(8) << transform[0] << endl;
        out << std::fixed << std::setprecision(8) << transform[3] << endl;
        out.close();

    }
    void point_trans(Eigen::Vector3d& in_pt, Eigen::Matrix4d T, double offx, double offy, Eigen::Vector3d& out_pt) {
        in_pt(0) -= offx;
        in_pt(1) -= offy;
        Eigen::Vector4d in_pt_homo(in_pt(0),in_pt(1),in_pt(2),1);
        Eigen::Vector4d out_pt_homo = T * in_pt_homo;
        out_pt(0) = out_pt_homo(0)  + offx;
        out_pt(1) = out_pt_homo(1)  + offy;
        out_pt(2) = out_pt_homo(2);
        in_pt(0) += offx;
        in_pt(1) += offy;
    }
};


void DSM_REG_v1(Parameters& par) {
    DSM_REG dsm_reg(par);
    dsm_reg.START("");
    //export transformation results
    std::string src_name = fs::path(par.src_path_).filename().string();
    std::string ref_name = fs::path(par.ref_path_).filename().string();
    src_name = src_name.substr(0,src_name.find(".tif"));
    ref_name = ref_name.substr(0,ref_name.find(".tif"));
    std::string results_name = src_name + "_" + ref_name + "_reg.txt";
    std::string results_path;
    if (par.out_dir_ == "") {
        results_path = fs::path(par.src_path_).replace_filename(results_name).string();
    }
    else {
        results_path = (fs::path(par.out_dir_) / results_name).string();
    }
    
    ofstream ofs;
    ofs.open(results_path);
    ofs << "### Final RMSE, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, T00, T01, T02, T03, T10, T11, T12, T13, T20, T21, T22, T23, T30, T31, T32, T33\n";
    ofs  << dsm_reg.final_rmse_ << "\n";
    ofs << std::setprecision(11) << dsm_reg.GLOBAL_OFFSET_X_m_ << "\n" << dsm_reg.GLOBAL_OFFSET_Y_m_ << "\n";
    ofs << std::setprecision(11) << dsm_reg.T_final_(0,0) << "\n" << dsm_reg.T_final_(0, 1) << "\n" << dsm_reg.T_final_(0, 2) << "\n" << dsm_reg.T_final_(0, 3) << "\n";
    ofs << std::setprecision(11) << dsm_reg.T_final_(1, 0) << "\n" << dsm_reg.T_final_(1, 1) << "\n" << dsm_reg.T_final_(1, 2) << "\n" << dsm_reg.T_final_(1, 3) << "\n";
    ofs << std::setprecision(11) << dsm_reg.T_final_(2, 0) << "\n" << dsm_reg.T_final_(2, 1) << "\n" << dsm_reg.T_final_(2, 2) << "\n" << dsm_reg.T_final_(2, 3) << "\n";
    ofs << std::setprecision(11) << dsm_reg.T_final_(3, 0) << "\n" << dsm_reg.T_final_(3, 1) << "\n" << dsm_reg.T_final_(3, 2) << "\n" << dsm_reg.T_final_(3, 3) << "\n";
    ofs.close();

    if (par.gen_data_) {
        //perform transform source dsm to reference
        std::string dsm_out_name = fs::path(par.src_path_).filename().string();
        dsm_out_name = dsm_out_name.replace(dsm_out_name.find(".tif"), sizeof(".tif") - 1, "_reg.tif");
        std::string dsm_out_path = fs::path(par.src_path_).replace_filename(dsm_out_name).string();
        DSM_TRANSFORM dsm_trans(par.src_path_, dsm_out_path, dsm_reg.T_final_, dsm_reg.GLOBAL_OFFSET_X_m_, dsm_reg.GLOBAL_OFFSET_Y_m_);
        dsm_trans.START();
    }

}

int main(int argc, char* argv[]) {
    argparse::ArgumentParser program("Large-scale DSM registration based on ICP, registered result will be in the same folder of src file");
    program.add_argument("-src").required().help("source/moving DSM file");
    program.add_argument("-dst").required().help("reference/fixed DSM file");
    program.add_argument("-outdir").default_value("").help("outdir where registration results stored");
    program.add_argument("-outlier_thresh").default_value(-1.0).help("outlier threshold for correspondences").scan<'g', double>();
    program.add_argument("-icp_num_pts").default_value(0.001).help("[default: 0.001] #pts used for ICP").scan<'g', double>();
    program.add_argument("-icp_max_iter").default_value(100).help("[default: 100] Rough registration using ICP: maximum number of iterations").scan<'d',int>();
    program.add_argument("-icp_rmse_threshold").default_value(1e-5).help("[default: 1e-5] RMSR threshold for early stop").scan<'g', double>();
    program.add_argument("-type").default_value("rigid").help("transformation type, can be (rigid, translation), 'rigid' contains 3DoF for rotation, 3FoF for translation, 'translation' only contains 3DoF ");
    program.add_argument("-gen_data").default_value(0).help("0: not transform, 1: transform source data").scan<'d', int>();

    try {
        program.parse_args(argc, argv);    // Example: ./main --color orange
    }
    catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    std::string src_file= program.get<std::string>("-src");
    std::string dst_file = program.get<std::string>("-dst");
    std::string out_dir = program.get<std::string>("-outdir");
    int icp_max_iter = program.get<int>("-icp_max_iter");
    int gen_data = program.get<int>("-gen_data");
    double icp_num_pts_ratio = program.get<double>("-icp_num_pts");
    double icp_rmse_threshold= program.get<double>("-icp_rmse_threshold");
    double outlier_thresh = program.get<double>("-outlier_thresh");
    std::string type = program.get<std::string>("-type");
    Parameters par;
    par.src_path_ = src_file;
    par.ref_path_ = dst_file;
    par.out_dir_ = out_dir;
    par.icp_num_pts_ratio_ = icp_num_pts_ratio;
    par.rough_icp_max_iter_ = icp_max_iter;
    par.rough_icp_rmse_threshold_ = icp_rmse_threshold;
    par.fine_icp_max_iter_ = icp_max_iter;
    par.fine_icp_rmse_threshold_ = icp_rmse_threshold;
    par.type_ = type;
    par.outlier_thresh_ = outlier_thresh;
    par.gen_data_ = gen_data;

    DSM_REG_v1(par);
}