/**
 * @file mesh.cpp
 * @brief Physical states and mesh storage.
 */
#include "mesh.h"
#include "table.h"
#include "utils.h"
#include "config.h"
#include "read_profile.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hd_core
{
namespace
{
    /**
     * @brief Validate the face coordinates of one mesh axis.
     * @param[in] faces All face coordinates, including ghost-region faces.
     * @param[in] count Expected number of cells along the axis.
     * @pre count is positive.
     * @throws std::invalid_argument If the face count is not count + 1, a coordinate is non-finite, or spacing is not finite and positive.
     */
    void validate_axis(const std::vector<double>& faces, int count)
    {
        if (faces.size() != static_cast<std::size_t>(count) + 1)
            throw std::invalid_argument("Face count must equal total cells + 1");

        for (std::size_t j = 0; j < faces.size(); ++j)
        {
            if (!std::isfinite(faces[j]))
                throw std::invalid_argument("Face coordinate must be finite");
            if (j > 0)
            {
                const double width = faces[j] - faces[j - 1];
                if (!std::isfinite(width) || width <= 0.0)
                    throw std::invalid_argument("Faces must have finite positive separation");
            }
        }
    }
    }

    int MeshBlock::ntotal(int d) const
    {
        if (d < 0 || d >= NDIM)
            throw std::out_of_range("Invalid mesh direction");

        if (nactive[d] <= 0 || nghost[d] < 0)
            throw std::invalid_argument("Invalid mesh size");
        if (nghost[d] > (std::numeric_limits<int>::max() - nactive[d]) / 2)
            throw std::overflow_error("Mesh size overflow");

        return nactive[d] + 2 * nghost[d];
    }
    Index MeshBlock::ntotal() const
    {
        Index result{};

        for (int d = 0; d < NDIM; ++d)
        {
            result[d] = ntotal(d);
        }

        return result;  
    }
    void MeshBlock::build_mesh(config::SimulationConfig& cfg)
    {
        nactive = cfg.nactive;
        nghost = cfg.nghost;
        geometry = parse_geometry(cfg.geometry);
        for(int d=0;d<hd_core::NDIM;d++)
        {
            bc_upper[d] = boundary::parse_bc(cfg.boundary_upper[d]);
            bc_lower[d] = boundary::parse_bc(cfg.boundary_lower[d]);
        }
        if(cfg.init_mesh_from_file)
        {
            config::Table table = config::read_table(cfg.init_mesh_file_name);
            std::array<std::string,3> axis; 
            switch (geometry)
            {
                case Geometry::Cartesian:
                    axis = {"x","y","z"};
                    break;
                case Geometry::Spherical:
                    axis = {"r","theta","phi"};
                    break;
                case Geometry::Cylindrical:
                    axis = {"r","theta","z"};
                    break;
                default:
                    break;
            }
            for(int d=0;d<NDIM;d++)
            {
                auto axis_name = axis[d];
                auto active_faces = table.column(axis_name);
                if(active_faces.size() - 1 != nactive[d])
                {
                    throw std::runtime_error("Mesh Initializing File Error: Axis" + std::to_string(d) + ", Number of faces incosistent with nactive in config file");
                }
                const int count = ntotal(d);
                validate_axis(active_faces, nactive[d]);  // monotonic increasing is checked.

                const int first = nghost[d];            // In cell it's the index of first active cell, and hence the lower face of first active cell.
                const int last = first + nactive[d];    // In cell it's the index of first right ghost cell, and for faces it's the lower face of that cell, hence the upper face of last active cell. 
                std::vector<double> coordinates(static_cast<std::size_t>(count) + 1);
                std::copy(active_faces.begin(), active_faces.end(), coordinates.begin() + first);
                //last - first + 1 faces are filled, i.e. nactive = 1

                const double dx_left  = coordinates[first + 1] - coordinates[first];
                const double dx_right = coordinates[last] - coordinates[last - 1];
                const bool reflect_left = bc_lower[d] == boundary::BoundaryCondition::Reflective;
                const bool reflect_right = bc_upper[d] == boundary::BoundaryCondition::Reflective;
                if ((reflect_left || reflect_right) && nghost[d] > nactive[d])
                    throw std::invalid_argument("Mirror boundary requires nghost <= nactive");

                for (int g = 1; g <= nghost[d]; g++)
                {
                    coordinates[first - g] = reflect_left
                        ? coordinates[first] - (coordinates[first + g] - coordinates[first])
                        : std::fma(-static_cast<double>(g), dx_left, coordinates[first]);
                    coordinates[last + g] = reflect_right
                        ? coordinates[last] + (coordinates[last] - coordinates[last - g])
                        : std::fma(static_cast<double>(g), dx_right, coordinates[last]);
                }
                //finally we got 2nghost + nactive + 1, which is ntotal +1
                validate_axis(coordinates, count);
                face[d].swap(coordinates);
                cfg.lower[d] = face[d][first];
                cfg.upper[d] = face[d][last];
            }
        }
        else
        {
            if(utils::iequals(cfg.init_mesh_uniform_type,"linear"))
            {
                for (int d = 0; d < NDIM; ++d)
                {
                    set_linear_uniform_faces(d, cfg.lower[d], cfg.upper[d]);
                }
            }
            else if(utils::iequals(cfg.init_mesh_uniform_type,"log"))
            {
                if(geometry != Geometry::Spherical && geometry != Geometry::Cylindrical)
                {
                    throw std::runtime_error("Log uniform can only work for radial axis of spherical or cylindrical coordinates");
                }
                set_log_uniform_faces(0,cfg.lower[0],cfg.upper[0]);

            }
            else
            {
                throw std::runtime_error("Wrong init_mesh_uniform_type, should be linear or log");
            }
        }
        allocate(); // Validates coordinates and geometry before allocating cells.
        build_geometry();
        if (cfg.init_profile_from_file)
        {
            config::read_profile(cfg.init_profile_file_name, *this);
        }
    }
    void MeshBlock::build_mesh(const config::SimulationConfig& cfg)
    {
        auto resolved = cfg;
        build_mesh(resolved);
    }
    void MeshBlock::build_geometry()
    {
        validate_geometry();

        // Volumes use the full cell layout, including ghosts.
        std::array<std::size_t, NDIM> cell_shape{};
        for (int d = 0; d < NDIM; ++d)
            cell_shape[d] = static_cast<std::size_t>(ntotal(d));
        //number of cells
        const auto storage_size = [](const auto& shape, std::size_t max_size) {
            std::size_t count = 1;
            for (const auto n : shape)
            {
                if (count > max_size / n)
                    throw std::length_error("Too many geometry entries");
                count *= n;
            }
            return count;
        };

        // Build temporary arrays so a failed calculation cannot partially replace
        // the old caches. No dependence on cells or on configuration bounds.
        std::vector<double> new_volumes;
        std::array<FaceField<double>, NDIM> new_areas;
        new_volumes.resize(storage_size(cell_shape, new_volumes.max_size()));
        const long double four_pi = 4.0L * std::acos(-1.0L);

        for (std::size_t flat = 0; flat < new_volumes.size(); ++flat)
        {
            long double volume = 1.0L;
            if (geometry == Geometry::Cartesian)
            {
                auto remainder = flat;
                for (int d = 0; d < NDIM; ++d)
                {
                    const auto i = remainder % cell_shape[d];
                    remainder /= cell_shape[d];
                    volume *= static_cast<long double>(face[d][i + 1]) - face[d][i];
                }
                //Use the similar method to get the position index, flat % x_length = index along i, then flat/= x_length, then flat % y_length to get y index, then ...
            }
            else // validate_geometry() restricts this to spherical 1D.
            {
                const long double lower = face[0][flat];
                const long double upper = face[0][flat + 1];
                // Factor the difference of cubes to avoid cancellation in thin
                // shells. Integrating r^2 also gives positive mirrored ghosts.
                volume = (four_pi / 3.0L) * (upper - lower)
                    * (lower * lower + lower * upper + upper * upper);
            }
            const double value = static_cast<double>(volume);
            if (!std::isfinite(value) || value <= 0.0)
                throw std::invalid_argument("Cell volume must be finite and positive");
            new_volumes[flat] = value;
        }

        for (int direction = 0; direction < NDIM; ++direction)
        {
            auto shape = ntotal();
            if (shape[direction] == std::numeric_limits<int>::max())
                throw std::overflow_error("Face shape exceeds the Index range");
            ++shape[direction];             //Shape of faces in current normal direction --> along normal direction there will be one more layer of faces
            auto& areas = new_areas[direction];     
            areas.allocate(shape); // Full face field, including ghost-region faces.
            for (std::size_t flat = 0; flat < areas.faces.size(); ++flat)
            {
                const auto index = areas.coordinates(flat);
                long double area = 1.0L;
                bool origin = false;
                if (geometry == Geometry::Cartesian)
                {
                    for (int d = 0; d < NDIM; ++d)
                    {
                        const auto i = index[d];
                        if (d != direction)             //Cartesian case, 1D area = 1.0, which is applicable when do A_LF_L-A_RF_R/V, which is equivalent to dF/dx. More than 1D, A_X = dydz, A_Y = ...(exlucd the direction) 
                            area *= static_cast<long double>(face[d][i + 1]) - face[d][i];
                    }
                }
                else
                {
                    const long double radius = face[0][index[0]];
                    area = four_pi * radius * radius;
                    origin = radius == 0.0L;
                }
                const double value = static_cast<double>(area);
                if (!std::isfinite(value) || value < 0.0 || (value == 0.0 && !origin))
                    throw std::invalid_argument("Face area must be finite and positive except at the origin");
                areas.faces[flat] = value;
            }
        }

        volume_cell.swap(new_volumes);
        face_area.swap(new_areas);
    }
    void MeshBlock::set_linear_uniform_faces(int d, double lower, double upper)
    {
        const int count = ntotal(d); // Checks direction, sizes, and integer overflow.
        if (!std::isfinite(lower) || !std::isfinite(upper) || upper <= lower)
            throw std::invalid_argument("Uniform grid bounds must be finite and increasing");

        const double length = upper - lower;
        const double width = length / nactive[d];
        if (!std::isfinite(length) || !std::isfinite(width) || width <= 0.0)
            throw std::invalid_argument("Invalid uniform grid spacing");

        std::vector<double> coordinates(static_cast<std::size_t>(count) + 1);
        const auto first = static_cast<std::size_t>(nghost[d]);     // In cell it's the index of first active cell, and hence the lower face of first active cell.
        const auto last = first + static_cast<std::size_t>(nactive[d]);     // In cell it's the index of first right ghost cell, and for faces it's the lower face of that cell, hence the upper face of last active cell. 
        for (std::size_t j = 0; j < coordinates.size(); ++j)
        {
            if (j == first)
                coordinates[j] = lower;
            else if (j == last)
                coordinates[j] = upper;
            else
                coordinates[j] = std::fma(static_cast<double>(j) - nghost[d], width, lower);
        }
        //nactive +1 +2nghost faces are set.
        // Detect overflow in ghost coordinates or collapsed faces due to rounding.
        validate_axis(coordinates, count);
        face[d].swap(coordinates); // Leave previous geometry unchanged on failure.
    }

    void MeshBlock::set_log_uniform_faces(int d, double lower, double upper)
    {
        const int count = ntotal(d); // Checks direction, sizes, and integer overflow.
        if (!std::isfinite(lower) || !std::isfinite(upper) || upper <= lower)
            throw std::invalid_argument("Uniform grid bounds must be finite and increasing");

        const double length = upper - lower;
        const double width = length / nactive[d];
        if (!std::isfinite(length) || !std::isfinite(width) || width <= 0.0)
            throw std::invalid_argument("Invalid uniform grid spacing");

        if(d!=0 || geometry == Geometry::Cartesian || lower<=0.0)
        {
            throw std::invalid_argument("Only radial faces in Spherical or Cylindrical Coordinates support log uniform");
        }
        std::vector<double> coordinates(static_cast<std::size_t>(count) + 1);
        const auto first = static_cast<std::size_t>(nghost[d]);
        const auto last = first + static_cast<std::size_t>(nactive[d]);
        for (std::size_t j = 0; j < coordinates.size(); ++j)
        {
            if (j == first)
                coordinates[j] = lower;
            else if (j == last)
                coordinates[j] = upper;
            else
            {
                const double t = (static_cast<double>(j) - nghost[d]) / nactive[d];
                coordinates[j] = lower * std::pow(upper/lower,t);
            }
            
        }
        // Detect overflow in ghost coordinates or collapsed faces due to rounding.
        validate_axis(coordinates, count);
        face[d].swap(coordinates); // Leave previous geometry unchanged on failure.
    }

    double MeshBlock::cell_width(int direction, int axis_index) const
    {
        // 只检查所选方向的单元数（含 ghost），不是整个网格的单元总数。
        const int axis_cell_count = ntotal(direction);
        if (axis_index < 0 || axis_index >= axis_cell_count)
            throw std::out_of_range("Cell index along the selected axis is out of range");
        if (face[direction].size() != static_cast<std::size_t>(axis_cell_count) + 1)
            throw std::invalid_argument("Face count must equal total cells + 1");

        // face[direction] 是该方向的一维面坐标数组，允许非均匀间距。
        // 沿该方向的第 axis_index 个单元由相邻两个面夹住。
        const auto lower_face_index = static_cast<std::size_t>(axis_index);
        const double lower_face_position = face[direction][lower_face_index];
        const double upper_face_position = face[direction][lower_face_index + 1];
        const double width = upper_face_position - lower_face_position;
        if (!std::isfinite(lower_face_position) || !std::isfinite(upper_face_position) ||
            !std::isfinite(width) || width <= 0.0)
            throw std::invalid_argument("Cell faces must have finite positive separation");
        return width;
    }

    double MeshBlock::cell_center(int direction, int axis_index) const
    {
        // cell_width() 同时检查方向、方向内下标及两侧面坐标是否合法。
        const double width = cell_width(direction, axis_index);
        const auto lower_face_index = static_cast<std::size_t>(axis_index);
        const double lower_face_position = face[direction][lower_face_index];
        // 返回所选方向上的坐标中点。
        return lower_face_position + 0.5 * width;
    }
    double MeshBlock::cell_volume_center(int direction, int axis_index) const
    {
        if(geometry==Geometry::Cartesian)
        {
            return cell_center(direction,axis_index);
        }
    if (geometry != Geometry::Spherical || NDIM != 1)
        throw std::invalid_argument(
            "cell_volume_center supports Cartesian and spherical 1D"
        );

        // validate
        (void)cell_width(direction, axis_index);
        const long double rL = face[direction][axis_index];
        const long double rR = face[direction][axis_index+1];   //though double, the intermediate result, e.g. r_L*r_L might be exceed the maximum of double,
                                                                // and after division the results might go back to double so here we use long double
        return static_cast<double> (0.75 * (rL+rR) * (rL*rL+rR*rR)/(rL*rL+rL*rR+rR*rR));
    }
    std::pair<double,double> MeshBlock::cell_face_pos(int direction, int axis_index) const
    {
        (void)cell_width(direction, axis_index);
        
        const double L = face[direction][axis_index];
        const double R = face[direction][axis_index+1];
        
        return {L,R};


    }
    double MeshBlock::cell_width(int direction, const Index& storage_index) const
    {
        if (direction < 0 || direction >= NDIM)
            throw std::out_of_range("Invalid mesh direction");
        // 张量积网格只需要该方向的坐标下标，无需转成 flat。
        return cell_width(direction, storage_index[direction]);
    }

    double MeshBlock::cell_center(int direction, const Index& storage_index) const
    {
        if (direction < 0 || direction >= NDIM)
            throw std::out_of_range("Invalid mesh direction");
        return cell_center(direction, storage_index[direction]);
    }

    double MeshBlock::cell_volume_center(int direction, const Index& storage_index) const
    {
        if (direction < 0 || direction >= NDIM)
            throw std::out_of_range("Invalid mesh direction");
        return cell_volume_center(direction, storage_index[direction]);
    }
    std::pair<double,double> MeshBlock::cell_face_pos(int direction, const Index& storage_index) const
    {
        if (direction < 0 || direction >= NDIM)
            throw std::out_of_range("Invalid mesh direction");
        return cell_face_pos(direction, storage_index[direction]);
    }
    
    double MeshBlock::active_length(int d) const
    {
        const int count = ntotal(d);
        validate_axis(face[d], count);
        const auto first = static_cast<std::size_t>(nghost[d]);
        const auto last = first + static_cast<std::size_t>(nactive[d]);
        const double length = face[d][last] - face[d][first];
        if (!std::isfinite(length) || length <= 0.0)
            throw std::invalid_argument("Active domain length must be finite and positive");
        return length;
    }

    void MeshBlock::validate_geometry() const
    {
        for (int d = 0; d < NDIM; ++d)
            (void)active_length(d); // Validates all faces, including ghost faces.

        // This is a pre-build coordinate check; caches may not exist yet.
        switch (geometry)
        {
            case Geometry::Cartesian:
                return;
            case Geometry::Spherical:
                if (NDIM != 1)
                    throw std::invalid_argument("Currently we only support 1D sphere");
                if (face[0][nghost[0]] < 0.0)
                    throw std::invalid_argument("Lower radial boundary can't be negative");
                return; // Negative ghost radii are allowed for mirrored states.
            case Geometry::Cylindrical:
                throw std::invalid_argument("Cylindrical geometry is not implemented");
        }
        throw std::invalid_argument("Invalid geometry");
    }

    void MeshBlock::validate_active() const
    {
        // Check dimensions and storage even when cells is empty.
        (void)flat_index(Index{}, true);

        for (std::size_t flat = 0; flat < cells.size(); ++flat)
        {
            const auto index = coordinates_active(flat);
            bool active = true;
            for (int d = 0; d < NDIM; ++d)
            {
                if (index[d] < 0 || index[d] >= nactive[d])
                {
                    active = false;
                    break;
                }
            }
            if (active)
                cells[flat].U.validate(cells[flat].gamma);
        }
    }

    void MeshBlock::allocate()
    {
        validate_geometry();
        std::size_t count = 1;

        for (int d = 0; d < NDIM; ++d)
        {
            const auto n = static_cast<std::size_t>(ntotal(d));

            if (count > cells.max_size() / n)
                throw std::length_error("Too many mesh cells");

            count *= n;
        }

        cells.resize(count);
    }

    Index MeshBlock::coordinates(std::size_t flat) const
    {
        if (flat >= cells.size())
            throw std::out_of_range("Cell index out of range");

        Index index{};
        std::size_t total = 1;

        for (int d = 0; d < NDIM; ++d)
        {
            const auto n = static_cast<std::size_t>(ntotal(d));
            if (total > cells.max_size() / n)
                throw std::length_error("Mesh is too large");

            total *= n;
            index[d] = static_cast<int>(flat % n);
            flat /= n;
        }

        if (cells.size() != total)
            throw std::logic_error("Cell storage does not match mesh dimensions");

        return index;
    }
    Index MeshBlock::coordinates_active(std::size_t flat) const
    {
        auto coordinates = MeshBlock::coordinates(flat);
        for(int i=0;i<NDIM;i++)
        {
            coordinates[i] -= nghost[i];
        }
        return coordinates;
    }

    std::size_t MeshBlock::flat_index(
        std::initializer_list<int> coordinates,
        bool active
    ) const
    {
        // 必须恰好提供 NDIM 个坐标。
        if (coordinates.size() != static_cast<std::size_t>(NDIM))
        {
            throw std::invalid_argument(
                "Number of coordinates must equal NDIM"
            );
        }

        Index index{};
        std::copy(coordinates.begin(), coordinates.end(), index.begin());
        return flat_index(index, active);
    }

    std::size_t MeshBlock::flat_index(
        const Index& coordinates,
        bool active
    ) const
    {
        std::size_t offset = 0;
        std::size_t stride = 1;
        int d = 0;

        for (int coordinate : coordinates)
        {
            const int n = ntotal(d);  // 包括 ghost cells
            const int limit = active ? nactive[d] : n;

            if (coordinate < 0 || coordinate >= limit)
            {
                throw std::out_of_range(
                    active ? "Active cell coordinate out of range"
                           : "Cell coordinate out of range"
                );
            }

            const auto count = static_cast<std::size_t>(n);

            // 防止总格数乘法溢出或超过容器容量限制。
            if (stride > cells.max_size() / count)
            {
                throw std::length_error("Mesh is too large");
            }

            // ntotal() has checked that the sizes and their sum fit in int.
            const int storage_coordinate =
                active ? coordinate + nghost[d] : coordinate;
            offset += static_cast<std::size_t>(storage_coordinate) * stride;
            stride *= count;
            ++d;
        }

        // 检查是否已按当前网格尺寸分配 cells。
        if (cells.size() != stride)
        {
            throw std::logic_error(
                "Cell storage does not match mesh dimensions"
            );
        }

        return offset;
    }

    Cell& MeshBlock::at(std::initializer_list<int> coordinates)
    {
        return cells.at(flat_index(coordinates));
    }

    const Cell& MeshBlock::at(
        std::initializer_list<int> coordinates
    ) const
    {
        return cells.at(flat_index(coordinates));
    }

    Cell& MeshBlock::at(const Index& coordinates)
    {
        return cells.at(flat_index(coordinates));
    }

    const Cell& MeshBlock::at(const Index& coordinates) const
    {
        return cells.at(flat_index(coordinates));
    }

    Cell& MeshBlock::at_active(std::initializer_list<int> coordinates)
    {
        return cells.at(flat_index(coordinates, true));
    }

    const Cell& MeshBlock::at_active(
        std::initializer_list<int> coordinates
    ) const
    {
        return cells.at(flat_index(coordinates, true));
    }

    Cell& MeshBlock::at_active(const Index& coordinates)
    {
        return cells.at(flat_index(coordinates, true));
    }

    const Cell& MeshBlock::at_active(const Index& coordinates) const
    {
        return cells.at(flat_index(coordinates, true));
    }

    Geometry parse_geometry(const std::string& name)
    {
        if(utils::iequals(name,"cartesian"))
        {
            return Geometry::Cartesian;
        }
        if(utils::iequals(name,"spherical"))
        {
            return Geometry::Spherical;
        }
        throw std::runtime_error("Invalid geometry");
    }
    bool MeshBlock::is_active(size_t flat) const
    {
        const auto index = coordinates_active(flat);
        bool active = true;

        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            if (index[d] < 0 || index[d] >= nactive[d])
            {
                active = false;
                break;
            }
        }
        return active;
    }
}
