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

using namespace std;
namespace fs = std::experimental::filesystem;
using namespace std::chrono;
#define PI 3.14159265;
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
    std::string ref_path_, src_path_;
    double valid_ratio_threshold_ = 0.5;
    double var_min_=4.0;
    int verbose_=1;
    // for General 
    double src_resolution_ = 0.5;
    double ref_resolution_ = 0.5;
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

    // for rough align
    int rough_num_pts_ = 24000;
    int search_half_size_pixel_ = 20;
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

void get_aoi_bbox_pixels(double* aoi_utm_bbox, double* ref_bbox, double* src_bbox, int ref_width, int ref_height, int src_width, int src_height, int* aoi_ref_bbox, int* aoi_src_bbox) {
    aoi_ref_bbox[0] = floor((aoi_utm_bbox[0] - ref_bbox[0]+0.25) / 0.5);
    aoi_ref_bbox[1] = floor((aoi_utm_bbox[1] - ref_bbox[0]+0.25) / 0.5);
    aoi_ref_bbox[2] = floor((ref_bbox[3]-aoi_utm_bbox[2]+ 0.25) / 0.5);
    aoi_ref_bbox[3] = floor((ref_bbox[3]-aoi_utm_bbox[3]+0.25) / 0.5);
    aoi_src_bbox[0] = floor((aoi_utm_bbox[0] - src_bbox[0]+0.25) / 0.5);
    aoi_src_bbox[1] = floor((aoi_utm_bbox[1] - src_bbox[0]+0.25) / 0.5);
    aoi_src_bbox[2] = floor((src_bbox[3]-aoi_utm_bbox[2]+0.25) / 0.5);
    aoi_src_bbox[3] = floor((src_bbox[3]-aoi_utm_bbox[3]+0.25) / 0.5);
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

void convert_ptcloud(double* data, int width, int height,double resolution,int resolution_factor, vector<gs::Point>& ptcloud, double offx, double offy) {
    ptcloud.clear();
    for (int row = 0; row < height; row+=resolution_factor) {
        for (int col = 0; col < width; col+=resolution_factor) {
            if (isnan(data[col + row * width])) {
                continue;
            }
            gs::Point pt((double)col * resolution+offx, (double)(- row) * resolution+offy, data[col + row * width]);
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
        poFeature->SetField("Ref_Var", cstr);
        delete[] cstr;

        var_str = to_string(tile_infos[i].src_var);
        char* cstr2 = new char[var_str.length() + 1];
        var_str.copy(cstr2, var_str.length());
        cstr2[var_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        poFeature->SetField("Src_Var", cstr2);
        delete[] cstr2;

        char* cstr1 = new char[id_str.length() + 1];
        id_str.copy(cstr1, id_str.length());
        cstr1[id_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        poFeature->SetField("Tile_ID", cstr1);
        delete[] cstr1;

        var_str = to_string(tile_infos[i].ref_valid);
        char* cstr3 = new char[var_str.length() + 1];
        var_str.copy(cstr3, var_str.length());
        cstr3[var_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        poFeature->SetField("Ref_Valid", cstr3);
        delete[] cstr3;

        var_str = to_string(tile_infos[i].src_valid);
        char* cstr4 = new char[var_str.length() + 1];
        var_str.copy(cstr4, var_str.length());
        cstr4[var_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        poFeature->SetField("Src_Valid", cstr4);
        delete[] cstr4;

        var_str = to_string(tile_infos[i].src_lambda1);
        char* cstr5 = new char[var_str.length() + 1];
        var_str.copy(cstr5, var_str.length());
        cstr5[var_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        poFeature->SetField("Src_Lmd1", cstr5);
        delete[] cstr5;

        var_str = to_string(tile_infos[i].src_lambda2);
        char* cstr6 = new char[var_str.length() + 1];
        var_str.copy(cstr6, var_str.length());
        cstr6[var_str.length()] = '\0';
        poFeature->SetField("Src_Lmd2", cstr6);
        delete[] cstr6;

        var_str = to_string(tile_infos[i].src_R);
        char* cstr7 = new char[var_str.length() + 1];
        var_str.copy(cstr7, var_str.length());
        cstr7[var_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        poFeature->SetField("Src_R", cstr7);
        delete[] cstr7;

        var_str = to_string(tile_infos[i].ref_lambda1);
        char* cstr8 = new char[var_str.length() + 1];
        var_str.copy(cstr8, var_str.length());
        cstr8[var_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        poFeature->SetField("Ref_Lmd1", cstr8);
        delete[] cstr8;

        var_str = to_string(tile_infos[i].ref_lambda2);
        char* cstr9 = new char[var_str.length() + 1];
        var_str.copy(cstr9, var_str.length());
        cstr9[var_str.length()] = '\0';
        poFeature->SetField("Ref_Lmd2", cstr9);
        delete[] cstr9;

        var_str = to_string(tile_infos[i].ref_R);
        char* cstr10 = new char[var_str.length() + 1];
        var_str.copy(cstr10, var_str.length());
        cstr10[var_str.length()] = '\0';
        //poFeature->SetField("Variance", cstr);
        poFeature->SetField("Ref_R", cstr10);
        delete[] cstr10;

        OGRGeometryH geo = OGR_G_CreateGeometry(wkbLinearRing);
        OGR_G_SetPoint_2D(geo, 0, tile_infos[i].utm_bbox[0], tile_infos[i].utm_bbox[2]);
        OGR_G_SetPoint_2D(geo, 1, tile_infos[i].utm_bbox[0], tile_infos[i].utm_bbox[3]);
        OGR_G_SetPoint_2D(geo, 2, tile_infos[i].utm_bbox[1], tile_infos[i].utm_bbox[3]);
        OGR_G_SetPoint_2D(geo, 3, tile_infos[i].utm_bbox[1], tile_infos[i].utm_bbox[2]);
        OGR_G_CloseRings(geo);
        geo = OGR_G_ForceToPolygon(geo);
        poFeature->SetGeometry(OGRGeometry::FromHandle(geo));
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

class DSM_TRANSFORM {
public:
    DSM_TRANSFORM(std::string in_path,std::string out_path,std::string ortho_path,std::string ortho_out_path,Eigen::Matrix4d Transform,double offx,double offy, double offz) {
        transform_ = Transform;
        in_zmin_ = -9999;
        in_zmax_ = -9999;
        int tile_size = 3000;

        double transform[6],transform_out[6];
        int xSize, ySize;
        out_path_ = out_path;
        in_path_ = in_path;
        GDALAllRegister();
        pdriver_ = GetGDALDriverManager()->GetDriverByName("GTiff");
        // dsm info
        in_dataset_ = (GDALDataset*)GDALOpen(in_path_.c_str(), GA_ReadOnly);
        in_dataset_->GetGeoTransform(transform);
        xSize = in_dataset_->GetRasterXSize();
        ySize = in_dataset_->GetRasterYSize();
        in_width_ = xSize;
        in_height_ = ySize;
        in_utm_bbox_[0] = transform[0]+transform[1]/2;
        in_utm_bbox_[1] = transform[0] + transform[1] / 2 + (xSize - 1) * transform[1];
        in_utm_bbox_[2] = transform[3]-transform[1] / 2 - (ySize - 1) * transform[1];
        in_utm_bbox_[3] = transform[3]- transform[1] / 2;

        GLOBAL_OFFSET_X_m_ = transform[0]+offx;
        GLOBAL_OFFSET_Y_m_ = transform[3]+offy;
        GLOBAL_OFFSET_Z_m_ = offz;

        resolution_ = transform[1];
        in_datatype_=GDALGetRasterDataType(GDALGetRasterBand(in_dataset_, 1));
        in_band_ = in_dataset_->GetRasterBand(1);

        // ortho info
        if (fs::exists(ortho_path)) {
            trans_ortho_ = true;
            double ortho_transform[6];
            int ortho_xSize, ortho_ySize;
            ortho_in_dataset_ = (GDALDataset*)GDALOpen(ortho_path.c_str(), GA_ReadOnly);
            ortho_in_dataset_->GetGeoTransform(ortho_transform);
            ortho_xSize = ortho_in_dataset_->GetRasterXSize();
            ortho_ySize = ortho_in_dataset_->GetRasterYSize();
            if (in_width_ != ortho_xSize || in_height_ != ortho_ySize || ortho_transform[1]!=resolution_) {
                trans_ortho_ = false;
                std::cout << "[WARNING] orthophoto height,width or gsd is not same as the dsm file, thus we only transform dsm file\n" << std::endl;
            }
        }
        else {
            trans_ortho_ = false;
        }

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
        point_trans(pt1, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, GLOBAL_OFFSET_Z_m_, pt1_trans);
        point_trans(pt2, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, GLOBAL_OFFSET_Z_m_, pt2_trans);
        point_trans(pt3, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, GLOBAL_OFFSET_Z_m_, pt3_trans);
        point_trans(pt4, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, GLOBAL_OFFSET_Z_m_, pt4_trans);
        point_trans(pt5, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, GLOBAL_OFFSET_Z_m_, pt5_trans);
        point_trans(pt6, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, GLOBAL_OFFSET_Z_m_, pt6_trans);
        point_trans(pt7, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, GLOBAL_OFFSET_Z_m_, pt7_trans);
        point_trans(pt8, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, GLOBAL_OFFSET_Z_m_, pt8_trans);
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
        //create empty output ortho
        if (trans_ortho_) {
            ortho_out_dataset_ = pdriver_->Create(ortho_out_path.c_str(), out_width_, out_height_, 3, GDT_Byte, NULL);
            ortho_out_dataset_->SetGeoTransform(transform_out);
            for (int band_id = 1; band_id < 4; ++band_id) {
                ortho_out_band_ = ortho_out_dataset_->GetRasterBand(band_id);
                num_tile_x = ceil((double)out_width_ / (double)tile_size);
                num_tile_y = ceil((double)out_height_ / (double)tile_size);
                for (int row = 0; row < num_tile_y; ++row) {
                    for (int col = 0; col < num_tile_x; ++col) {
                        fill_x = (tile_size * (col + 1) <= out_width_) ? tile_size * (col + 1) : out_width_;
                        fill_x -= tile_size * col;
                        fill_y = (tile_size * (row + 1) <= out_height_) ? tile_size * (row + 1) : out_height_;
                        fill_y -= tile_size * row;
                        rasterio_write_nan(ortho_out_band_, tile_size * col, tile_size * row, fill_x, fill_y);
                    }
                }
            }

            std::string out_tfw_path = ortho_out_path;
            strutil::replace_last(out_tfw_path, ".tif", ".tfw");
            export_tfw(out_tfw_path, transform_out);

        }
        //GDALClose(out_dataset_);
    };
    std::string in_path_, out_path_;
    GDALDriver* pdriver_;
    GDALDataset* out_dataset_;
    GDALDataset* in_dataset_;
    GDALDataset* ortho_out_dataset_;
    GDALDataset* ortho_in_dataset_;
    GDALRasterBand* in_band_;
    GDALRasterBand* out_band_;
    GDALRasterBand* ortho_in_band_;
    GDALRasterBand* ortho_out_band_;
    GDALDataType in_datatype_;
    int out_width_, out_height_, in_width_, in_height_;
    double GLOBAL_OFFSET_X_m_;
    double GLOBAL_OFFSET_Y_m_;
    double GLOBAL_OFFSET_Z_m_;
    double resolution_;
    double out_utm_bbox_[4];
    double in_utm_bbox_[4];
    double in_zmin_, in_zmax_;
    bool trans_ortho_ = true;
    Eigen::Matrix4d transform_;

    void START() {
        double* write_value = new double[1];
        double* read_value = new double[1];
        unsigned char* rgb = new unsigned char[3];
        //double pt[3],pt_out[3];
        Eigen::Vector3d pt, pt_out;
        int col_out, row_out;
        for (int row = 0; row < in_height_; ++row)
        {
            for (int col = 0; col < in_width_; ++col) {
                pt[0] = in_utm_bbox_[0] + resolution_ * col;
                pt[1] = in_utm_bbox_[3] - resolution_ * row;
                in_band_->RasterIO(GF_Read, col, row, 1, 1,
                    read_value, 1, 1, GDT_Float64, 0, 0);
                if (trans_ortho_) {
                    ortho_in_band_ = ortho_in_dataset_->GetRasterBand(1);
                    ortho_in_band_->RasterIO(GF_Read, col, row, 1, 1,&rgb[0], 1, 1, GDT_Byte, 0, 0);
                    ortho_in_band_ = ortho_in_dataset_->GetRasterBand(2);
                    ortho_in_band_->RasterIO(GF_Read, col, row, 1, 1, &rgb[1], 1, 1, GDT_Byte, 0, 0);
                    ortho_in_band_ = ortho_in_dataset_->GetRasterBand(3);
                    ortho_in_band_->RasterIO(GF_Read, col, row, 1, 1, &rgb[2], 1, 1, GDT_Byte, 0, 0);
                }
                pt[2] = read_value[0];
                if (isnan(pt[2])) {
                    continue;
                }
                point_trans(pt, transform_, GLOBAL_OFFSET_X_m_, GLOBAL_OFFSET_Y_m_, GLOBAL_OFFSET_Z_m_, pt_out);
                col_out = floor((pt_out[0] - out_utm_bbox_[0] + resolution_ / 2) / resolution_);
                row_out= floor((out_utm_bbox_[3]- pt_out[1]+ resolution_ / 2) / resolution_);
                if (row_out < 0) {
                    int a = 0;
                }
                //std::cout << col_out << "," << row_out << std::endl;
                out_band_->RasterIO(GF_Read, col_out, row_out, 1, 1,
                    read_value, 1, 1, GDT_Float64, 0, 0);
                if (pt_out[2] > read_value[0] || isnan(read_value[0])) {
                    out_band_->RasterIO(GF_Write, col_out, row_out, 1, 1,
                        &pt_out[2], 1, 1, GDT_Float64, 0, 0);
                    if (trans_ortho_) {
                        ortho_out_band_ = ortho_out_dataset_->GetRasterBand(1);
                        ortho_out_band_->RasterIO(GF_Write, col_out, row_out, 1, 1, &rgb[0], 1, 1, GDT_Byte, 0, 0);
                        ortho_out_band_ = ortho_out_dataset_->GetRasterBand(2);
                        ortho_out_band_->RasterIO(GF_Write, col_out, row_out, 1, 1, &rgb[1], 1, 1, GDT_Byte, 0, 0);
                        ortho_out_band_ = ortho_out_dataset_->GetRasterBand(3);
                        ortho_out_band_->RasterIO(GF_Write, col_out, row_out, 1, 1, &rgb[2], 1, 1, GDT_Byte, 0, 0);
                    }
                }
            }
        }
        GDALClose(in_dataset_);
        GDALClose(out_dataset_);
        if (trans_ortho_) {
            GDALClose(ortho_in_dataset_);
            GDALClose(ortho_out_dataset_);
        }
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
        case 8: {
            //GInt16* buf16 = new GInt16[cnt];
            unsigned char* buf8 = (unsigned char*)CPLMalloc(sizeof(unsigned char) * cnt);
            for (int i = 0; i < cnt; ++i) {
                buf8[i] = std::numeric_limits<unsigned char>::quiet_NaN();
            }
            band->RasterIO(GF_Write, offsetx, offsety, xsize, ysize, buf8, xsize, ysize, GDT_Byte, 0, 0);
            //delete[] buf16;
            CPLFree(buf8);
            break;
        }
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
    void point_trans(Eigen::Vector3d& in_pt, Eigen::Matrix4d T, double offx, double offy, double offz, Eigen::Vector3d& out_pt) {
        in_pt(0) -= offx;
        in_pt(1) -= offy;
        in_pt(2) -= offz;
        Eigen::Vector4d in_pt_homo(in_pt(0),in_pt(1),in_pt(2),1);
        Eigen::Vector4d out_pt_homo = T * in_pt_homo;
        out_pt(0) = out_pt_homo(0)  + offx;
        out_pt(1) = out_pt_homo(1)  + offy;
        out_pt(2) = out_pt_homo(2) + offz;
        in_pt(0) += offx;
        in_pt(1) += offy;
        in_pt(2) += offz;
    }
};


int main(int argc, char* argv[]) {
    argparse::ArgumentParser program("Apply transform to DSM, transformed result will be in the same folder of input file, ");
    program.add_argument("-dsm").required().help("input DSM file, should come with tfw file");
    program.add_argument("-ortho").default_value("").help(" (optional) input orthophoto file, should be in the same gsd,height,width as dsm file, and come with tfw file");
    program.add_argument("-outdir").default_value("").help(" (optional) outdir where transformed dsm stored");
    program.add_argument("-rotation").required().help("rotation parameters 9 elements, r11,r12,r13,r21,r22,r23,r31,r32,r33").nargs(9).scan<'g',double>();
    program.add_argument("-translation").required().help("translation parameters 3 elements, t1,t2,t3").nargs(3).scan<'g', double>();

    program.add_argument("-reso").default_value(0.5).help("DSM resolution").scan<'g', double>();
    program.add_argument("-rot_center_x").default_value(0.0).help("Note: in local coordinate system, rotation center in x direction, 0 represents lefttop corner of DSM").scan<'g', double>();
    program.add_argument("-rot_center_y").default_value(0.0).help("Note: in local coordinate system, rotation center in y direction, 0 represents lefttop corner of DSM").scan<'g', double>();
    program.add_argument("-rot_center_z").default_value(0.0).help("Note: in local coordinate system, rotation center in z direction, 0 represents lefttop corner of DSM").scan<'g', double>();

    try {
        program.parse_args(argc, argv);    // Example: ./main --color orange
    }
    catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    std::string input_file= program.get<std::string>("-dsm");
    std::string ortho_file = program.get<std::string>("-ortho");
    std::string outdir = program.get<std::string>("-outdir");

    double reso = program.get<double>("-reso");
    double r_offx = program.get<double>("-rot_center_x");
    double r_offy = program.get<double>("-rot_center_y");
    double r_offz = program.get<double>("-rot_center_z");
    std::vector<double> r = program.get<std::vector<double>>("-rotation");
    std::vector<double> t = program.get<std::vector<double>>("-translation");
    Eigen::Matrix4d T;
    T << r[0], r[1], r[2], t[0],
        r[3], r[4], r[5], t[1],
        r[6], r[7], r[8], t[2],
        0, 0, 0, 1;

    std::string dsm_out = input_file;
    if (outdir == "") {
        dsm_out = dsm_out.replace(dsm_out.find(".tif"), sizeof(".tif") - 1, "_transformed.tif");
    }
    else {
        std::string dsm_out_name = fs::path(input_file).filename().string();
        dsm_out_name= dsm_out_name.replace(dsm_out_name.find(".tif"), sizeof(".tif") - 1, "_transformed.tif");
        dsm_out = (fs::path(outdir)/dsm_out_name).string();
    }


    std::string ortho_out = ortho_file;
    if (ortho_out != "") {
        if (outdir == "") {
            ortho_out = ortho_out.replace(ortho_out.find(".tif"), sizeof(".tif") - 1, "_transformed.tif");
        }
        else {
            std::string ortho_out_name = fs::path(ortho_file).filename().string();
            ortho_out_name = ortho_out_name.replace(ortho_out_name.find(".tif"), sizeof(".tif") - 1, "_transformed.tif");
            ortho_out = fs::path(outdir).replace_filename(ortho_out_name).string();
        }
        
    }


    DSM_TRANSFORM dsm_trans(input_file,dsm_out,ortho_file,ortho_out, T, r_offx, r_offy, r_offz);
    dsm_trans.START();
	
}
