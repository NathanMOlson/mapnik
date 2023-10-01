/*****************************************************************************
 *
 * This file is part of Mapnik (c++ mapping toolkit)
 *
 * Copyright (C) 2021 Artem Pavlenko
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *****************************************************************************/

// mapnik
#include <mapnik/geometry/boost_adapters.hpp>
#include <mapnik/geometry/box2d.hpp>
#include <mapnik/geometry/multi_point.hpp>
#include <mapnik/projection.hpp>
#include <mapnik/proj_transform.hpp>
#include <mapnik/coord.hpp>
#include <mapnik/util/is_clockwise.hpp>
#include <mapnik/util/trim.hpp>
// boost
#include <boost/geometry/algorithms/envelope.hpp>

#ifdef MAPNIK_USE_PROJ
// proj
#include <proj.h>
#endif

// stl
#include <vector>
#include <stdexcept>

namespace mapnik {

namespace { // (local)

// Returns points in clockwise order. This allows us to do anti-meridian checks.
template<typename T>
auto envelope_points(box2d<T> const& env, std::size_t num_points) -> geometry::multi_point<T>
{
    auto width = env.width();
    auto height = env.height();

    geometry::multi_point<T> coords;
    coords.reserve(num_points);

    // top side: left >>> right
    // gets extra point if (num_points % 4 >= 1)
    for (std::size_t i = 0, n = (num_points + 3) / 4; i < n; ++i)
    {
        auto x = env.minx() + (i * width) / n;
        coords.emplace_back(x, env.maxy());
    }

    // right side: top >>> bottom
    // gets extra point if (num_points % 4 >= 3)
    for (std::size_t i = 0, n = (num_points + 1) / 4; i < n; ++i)
    {
        auto y = env.maxy() - (i * height) / n;
        coords.emplace_back(env.maxx(), y);
    }

    // bottom side: right >>> left
    // gets extra point if (num_points % 4 >= 2)
    for (std::size_t i = 0, n = (num_points + 2) / 4; i < n; ++i)
    {
        auto x = env.maxx() - (i * width) / n;
        coords.emplace_back(x, env.miny());
    }

    // left side: bottom >>> top
    // never gets extra point
    for (std::size_t i = 0, n = (num_points + 0) / 4; i < n; ++i)
    {
        auto y = env.miny() + (i * height) / n;
        coords.emplace_back(env.minx(), y);
    }

    return coords;
}

} // namespace

proj_transform::proj_transform(projection const& source, projection const& dest)
    : is_source_longlat_(false)
    , is_dest_longlat_(false)
    , is_source_equal_dest_(false)
    , wgs84_to_merc_(false)
    , merc_to_wgs84_(false)
    , is_dest_camera_(false)
    , is_src_camera_(false)
{
    is_source_equal_dest_ = (source == dest);
    if (!is_source_equal_dest_)
    {
        is_source_longlat_ = source.is_geographic();
        is_dest_longlat_ = dest.is_geographic();
        boost::optional<well_known_srs_e> src_k = source.well_known();
        boost::optional<well_known_srs_e> dest_k = dest.well_known();
        bool known_trans = false;
        projection src_proj = source;
        projection dest_proj = dest;
        if (source.is_camera())
        {
            is_src_camera_ = true;
            init_cam_params(source.params());
            if(dest_k && *dest_k == well_known_srs_enum::WGS_84)
            {
                known_trans = true;
            }
            else
            {
                src_proj = projection("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs");
            }
        }
        if (dest.is_camera())
        {
            is_dest_camera_ = true;
            init_cam_params(dest.params());
            if(src_k && *src_k == well_known_srs_enum::WGS_84)
            {
                known_trans = true;
            }
            else
            {
                dest_proj = projection("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs");
            }
        }
        if (src_k && dest_k)
        {
            if (*src_k == well_known_srs_enum::WGS_84 && *dest_k == well_known_srs_enum::WEB_MERC)
            {
                wgs84_to_merc_ = true;
                known_trans = true;
            }
            else if (*src_k == well_known_srs_enum::WEB_MERC && *dest_k == well_known_srs_enum::WGS_84)
            {
                merc_to_wgs84_ = true;
                known_trans = true;
            }
        }
        if (!known_trans)
        {
#ifdef MAPNIK_USE_PROJ
            ctx_ = proj_context_create();
            proj_log_level(ctx_, PJ_LOG_ERROR);
            transform_ = proj_create_crs_to_crs(ctx_, src_proj.params().c_str(), dest_proj.params().c_str(), nullptr);
            if (transform_ == nullptr)
            {
                throw std::runtime_error(
                  std::string("Cannot initialize proj_transform  (crs_to_crs) for given projections: '") +
                  src_proj.params() + "'->'" + dest_proj.params() +
#if MAPNIK_PROJ_VERSION >= 80000
                  "' because of " + std::string(proj_context_errno_string(ctx_, proj_context_errno(ctx_))));
#else
                  "'");
#endif
            }
            PJ* transform_gis = proj_normalize_for_visualization(ctx_, transform_);
            if (transform_gis == nullptr)
            {
                throw std::runtime_error(
                  std::string("Cannot initialize proj_transform (normalize) for given projections: '") +
                  src_proj.params() + "'->'" + dest_proj.params() +
#if MAPNIK_PROJ_VERSION >= 80000
                  "' because of " + std::string(proj_context_errno_string(ctx_, proj_context_errno(ctx_))));
#else
                  "'");
#endif
            }
            proj_destroy(transform_);
            transform_ = transform_gis;
#else
            throw std::runtime_error(
              std::string(
                "Cannot initialize proj_transform for given projections without proj support (-DMAPNIK_USE_PROJ): '") +
              src_proj.params() + "'->'" + dest_proj.params() + "'");
#endif
        }
    }
}

proj_transform::~proj_transform()
{
#ifdef MAPNIK_USE_PROJ
    if (transform_)
    {
        proj_destroy(transform_);
        transform_ = nullptr;
    }
    if (ctx_)
    {
        proj_context_destroy(ctx_);
        ctx_ = nullptr;
    }
#endif
}

bool proj_transform::equal() const
{
    return is_source_equal_dest_;
}

bool proj_transform::is_known() const
{
    return merc_to_wgs84_ || wgs84_to_merc_;
}

bool proj_transform::forward(double& x, double& y, double& z) const
{
    return forward(&x, &y, &z, 1);
}

bool proj_transform::forward(geometry::point<double>& p) const
{
    double z = 0;
    return forward(&(p.x), &(p.y), &z, 1);
}

unsigned int proj_transform::forward(std::vector<geometry::point<double>>& ls) const
{
    std::size_t size = ls.size();
    if (size == 0)
        return 0;

    if (is_source_equal_dest_)
        return 0;

    if (wgs84_to_merc_)
    {
        lonlat2merc(ls);
        return 0;
    }
    else if (merc_to_wgs84_)
    {
        merc2lonlat(ls);
        return 0;
    }

    geometry::point<double>* ptr = ls.data();
    double* x = reinterpret_cast<double*>(ptr);
    double* y = x + 1;
    double* z = nullptr;
    if (!forward(x, y, z, size, 2))
    {
        return size;
    }
    return 0;
}

bool proj_transform::forward(double* x, double* y, double* z, std::size_t point_count, std::size_t offset) const
{
    if (is_source_equal_dest_)
        return true;

    if (wgs84_to_merc_)
    {
        return lonlat2merc(x, y, point_count, offset);
    }
    else if (merc_to_wgs84_)
    {
        return merc2lonlat(x, y, point_count, offset);
    }

#ifdef MAPNIK_USE_PROJ
    if (transform_)
    {
    if (proj_trans_generic(transform_,
                           PJ_FWD,
                           x,
                           offset * sizeof(double),
                           point_count,
                           y,
                           offset * sizeof(double),
                           point_count,
                           z,
                           offset * sizeof(double),
                           point_count,
                           0,
                           0,
                           0) != point_count)
        return false;
    }

#endif
    if (is_dest_camera_)
    {
        return lonlat2camera(x, y, z, point_count, offset);
    }
    if (is_src_camera_)
    {
        return camera2latlon(x, y, z, point_count, offset);
    }
    return true;
}

bool proj_transform::backward(double* x, double* y, double* z, std::size_t point_count, std::size_t offset) const
{
    if (is_source_equal_dest_)
        return true;

    if (wgs84_to_merc_)
    {
        return merc2lonlat(x, y, point_count, offset);
    }
    else if (merc_to_wgs84_)
    {
        return lonlat2merc(x, y, point_count, offset);
    }

#ifdef MAPNIK_USE_PROJ
    if (transform_)
    {
    if (proj_trans_generic(transform_,
                           PJ_INV,
                           x,
                           offset * sizeof(double),
                           point_count,
                           y,
                           offset * sizeof(double),
                           point_count,
                           z,
                           offset * sizeof(double),
                           point_count,
                           0,
                           0,
                           0) != point_count)
        return false;
    }
#endif
    if (is_dest_camera_)
    {
        return camera2latlon(x, y, z, point_count, offset);
    }
    if (is_src_camera_)
    {
        return lonlat2camera(x, y, z, point_count, offset);
    }
    return true;
}

bool proj_transform::backward(double& x, double& y, double& z) const
{
    return backward(&x, &y, &z, 1);
}

bool proj_transform::backward(geometry::point<double>& p) const
{
    double z = 0;
    return backward(&(p.x), &(p.y), &z, 1);
}

unsigned int proj_transform::backward(std::vector<geometry::point<double>>& ls) const
{
    std::size_t size = ls.size();
    if (size == 0)
        return 0;

    if (is_source_equal_dest_)
        return 0;

    if (wgs84_to_merc_)
    {
        merc2lonlat(ls);
        return 0;
    }
    else if (merc_to_wgs84_)
    {
        lonlat2merc(ls);
        return 0;
    }

    geometry::point<double>* ptr = ls.data();
    double* x = reinterpret_cast<double*>(ptr);
    double* y = x + 1;
    double* z = nullptr;
    if (!backward(x, y, z, size, 2))
    {
        return size;
    }
    return 0;
}

bool proj_transform::forward(box2d<double>& box) const
{
    if (is_source_equal_dest_)
        return true;

    double llx = box.minx();
    double ulx = box.minx();
    double lly = box.miny();
    double lry = box.miny();
    double lrx = box.maxx();
    double urx = box.maxx();
    double uly = box.maxy();
    double ury = box.maxy();
    double z = 0.0;
    if (!forward(llx, lly, z))
        return false;
    if (!forward(lrx, lry, z))
        return false;
    if (!forward(ulx, uly, z))
        return false;
    if (!forward(urx, ury, z))
        return false;

    double minx = std::min(ulx, llx);
    double miny = std::min(lly, lry);
    double maxx = std::max(urx, lrx);
    double maxy = std::max(ury, uly);
    box.init(minx, miny, maxx, maxy);
    return true;
}

bool proj_transform::backward(box2d<double>& box) const
{
    if (is_source_equal_dest_)
        return true;

    double x[4], y[4];
    x[0] = box.minx(); // llx 0
    y[0] = box.miny(); // lly 1
    x[1] = box.maxx(); // lrx 2
    y[1] = box.miny(); // lry 3
    x[2] = box.minx(); // ulx 4
    y[2] = box.maxy(); // uly 5
    x[3] = box.maxx(); // urx 6
    y[3] = box.maxy(); // ury 7

    if (!backward(x, y, nullptr, 4, 1))
        return false;

    double minx = std::min(x[0], x[2]);
    double miny = std::min(y[0], y[1]);
    double maxx = std::max(x[1], x[3]);
    double maxy = std::max(y[2], y[3]);
    box.init(minx, miny, maxx, maxy);
    return true;
}

// More robust, but expensive, bbox transform
// in the face of proj4 out of bounds conditions.
// Can result in 20 -> 10 r/s performance hit.
// Alternative is to provide proper clipping box
// in the target srs by setting map 'maximum-extent'

bool proj_transform::backward(box2d<double>& env, std::size_t points) const
{
    if (is_source_equal_dest_)
        return true;

    if (wgs84_to_merc_ || merc_to_wgs84_)
    {
        return backward(env);
    }

    auto coords = envelope_points(env, points); // this is always clockwise

    for (auto& p : coords)
    {
        double z = 0;
        if (!backward(p.x, p.y, z))
            return false;
    }

    box2d<double> result;
    boost::geometry::envelope(coords, result);

    if (is_source_longlat_ && !util::is_clockwise(coords))
    {
        // we've gone to a geographic CS, and our clockwise envelope has
        // changed into an anticlockwise one. This means we've crossed the antimeridian, and
        // need to expand the X direction to +/-180 to include all the data. Once we can deal
        // with multiple bboxes in queries we can improve.
        double miny = result.miny();
        result.expand_to_include(-180.0, miny);
        result.expand_to_include(180.0, miny);
    }

    env.re_center(result.center().x, result.center().y);
    env.height(result.height());
    env.width(result.width());
    return true;
}

bool proj_transform::forward(box2d<double>& env, std::size_t points) const
{
    if (is_source_equal_dest_)
        return true;

    if (wgs84_to_merc_ || merc_to_wgs84_)
    {
        return forward(env);
    }

    auto coords = envelope_points(env, points); // this is always clockwise

    for (auto& p : coords)
    {
        double z = 0;
        if (!forward(p.x, p.y, z))
            return false;
    }

    box2d<double> result;
    boost::geometry::envelope(coords, result);

    if (is_dest_longlat_ && !util::is_clockwise(coords))
    {
        // we've gone to a geographic CS, and our clockwise envelope has
        // changed into an anticlockwise one. This means we've crossed the antimeridian, and
        // need to expand the X direction to +/-180 to include all the data. Once we can deal
        // with multiple bboxes in queries we can improve.
        double miny = result.miny();
        result.expand_to_include(-180.0, miny);
        result.expand_to_include(180.0, miny);
    }

    env.re_center(result.center().x, result.center().y);
    env.height(result.height());
    env.width(result.width());

    return true;
}

std::string proj_transform::definition() const
{
#ifdef MAPNIK_USE_PROJ
    if (transform_)
    {
        PJ_PROJ_INFO info = proj_pj_info(transform_);
        return mapnik::util::trim_copy(info.definition);
    }
    else
#endif
      if (wgs84_to_merc_)
    {
        return "wgs84 => merc";
    }
    else if (merc_to_wgs84_)
    {
        return "merc => wgs84";
    }
    return "unknown";
}

void proj_transform::init_cam_params(const std::string& params)
{
    double roll = 0;
    double pitch = 0;
    double yaw = 0;

    std::stringstream ss;
    
    int index = params.find("+lat=") + 5;
    ss = std::stringstream(params.substr(index, std::string::npos));
    ss >> cam_params_.lat;
    
    index = params.find("+lon=") + 5;
    ss = std::stringstream(params.substr(index, std::string::npos));
    ss >> cam_params_.lon;
    
    index = params.find("+alt=") + 5;
    ss = std::stringstream(params.substr(index, std::string::npos));
    ss >> cam_params_.alt;
    
    index = params.find("+ifov=") + 6;
    ss = std::stringstream(params.substr(index, std::string::npos));
    ss >> cam_params_.i_fov;
    
    index = params.find("+width=") + 7;
    ss = std::stringstream(params.substr(index, std::string::npos));
    ss >> cam_params_.width;

    index = params.find("+height=") + 8;
    ss = std::stringstream(params.substr(index, std::string::npos));
    ss >> cam_params_.height;

    

    index = params.find("+roll=");
    if(index != std::string::npos)
    {
        index += 6;
        ss = std::stringstream(params.substr(index, std::string::npos));
        ss >> roll;
    }

    index = params.find("+pitch=");
    if(index != std::string::npos)
    {
        index += 7;
        ss = std::stringstream(params.substr(index, std::string::npos));
        ss >> pitch;
    }

    index = params.find("+yaw=");
    if(index != std::string::npos)
    {
        index += 5;
        ss = std::stringstream(params.substr(index, std::string::npos));
        ss >> yaw;
    }

    double croll = cos(util::radians(roll));
    double sroll = sin(util::radians(roll));
    double cpitch = cos(util::radians(pitch));
    double spitch = sin(util::radians(pitch));
    double cyaw = cos(util::radians(yaw));
    double syaw = sin(util::radians(yaw));

    cam_params_.rot[0][0] = cyaw*cpitch;
    cam_params_.rot[0][1] = syaw*cpitch;
    cam_params_.rot[0][2] = -spitch;
    cam_params_.rot[1][0] = cyaw*spitch*sroll - syaw*croll;
    cam_params_.rot[1][1] = syaw*spitch*sroll + cyaw*croll;
    cam_params_.rot[1][2] = cpitch*sroll;
    cam_params_.rot[2][0] = cyaw*spitch*croll + syaw*sroll;
    cam_params_.rot[2][1] = syaw*spitch*croll - cyaw*sroll;
    cam_params_.rot[2][2] = cpitch*croll;
}

bool proj_transform::lonlat2camera(double* x, double* y, const double* z, std::size_t point_count, std::size_t stride) const
{
    bool valid = true;
    for (std::size_t i = 0; i < point_count; ++i)
    {
        double e = EARTH_RADIUS * util::radians(x[i * stride] - cam_params_.lon) * cos(util::radians(cam_params_.lat));
        double n = EARTH_RADIUS * util::radians(y[i * stride] - cam_params_.lat);
        double d = cam_params_.alt - z[i * stride];
        double f = n*cam_params_.rot[0][0] + e*cam_params_.rot[0][1] + d*cam_params_.rot[0][2];
        double r = n*cam_params_.rot[1][0] + e*cam_params_.rot[1][1] + d*cam_params_.rot[1][2];
        double b = n*cam_params_.rot[2][0] + e*cam_params_.rot[2][1] + d*cam_params_.rot[2][2];
        x[i * stride] = (cam_params_.width - 1)/2.0 + r / b / cam_params_.i_fov;
        y[i * stride] = (cam_params_.height - 1)/2.0 + f / b / cam_params_.i_fov;
        if(b < 0)
        {
            valid = false;
        }
    }
    return valid;
}

bool proj_transform::camera2latlon(double* x, double* y, double* z, std::size_t point_count, std::size_t stride) const
{
    for (std::size_t i = 0; i < point_count; ++i)
    {
        double d = cam_params_.alt;
        double e = d * cam_params_.i_fov * (x[i * stride] - (cam_params_.width - 1)/2.0);
        double n = d * cam_params_.i_fov * (y[i * stride] - (cam_params_.height - 1)/2.0);
        x[i * stride] = cam_params_.lon + util::degrees(e / EARTH_RADIUS) / cos(util::radians(cam_params_.lat));
        y[i * stride] = cam_params_.lat + util::degrees(n / EARTH_RADIUS);
        z[i * stride] = 0;
    }
    return true;
}

} // namespace mapnik
