#pragma once
#include "elements_io.hpp"

static inline void mesh_settings() {
  cout << "C - Cell size:            " << mesh.cell_size << endl
       << "X - Origin x:             " << mesh.x0 << endl
       << "Y - Origin y:             " << mesh.y0 << endl
       << "Z - Origin z:             " << mesh.z0 << endl;
  char opt = 'd';
  cin >> opt;
  switch (opt) {
  case 'C': {
    cin >> mesh.cell_size;
  } break;
  case 'X': {
    cin >> mesh.x0;
  } break;
  case 'Y': {
    cin >> mesh.y0;
  } break;
  case 'Z': {
    cin >> mesh.z0;
  } break;
  default:
    COMMAND_NOT_FOUND;
    break;
  }
}

template <typename Callable>
static inline void input_patch_indices(Callable &&callable, u32 max) {
  while (true) {
    u32 i = max;
    cin >> i;
    if (i < max)
      callable(i);
  }
}

static inline void select_patches() {
  while (true) {
    cout << "T - Select some patches.\n"
            "F - Unselect some patches.\n"
            "R - Reverse selection of some patches.\n"
            "d - discard.\n";
    char opt = 'd';
    switch (opt) {
    case 'T':
      break;
    case 'F':
      break;
    case 'R':
      break;
    case 'd':
      return;
    default:
      COMMAND_NOT_FOUND;
      break;
    }
  }
}

static inline void attach(const vector<patch_info> &patches,
                          const vector<frame> &frames,
                          const tuple<vector<float>, vector<u32>, vector<u32>,
                                      vector<u32>, vector<u32>> &element_data) {
  auto &[nodes, element_sizes, element_indices, node_numbers, element_numbers] =
      element_data;

  auto elem_avail = [&]() {
    auto &[nodes, element_sizes, element_indices, node_numbers,
           element_numbers] = element_data;

    return !nodes.empty() && !element_sizes.empty() && !element_indices.empty();
  };

  /// polygon sizes, polygon indices, element sizes, node number, element
  /// number, surface numbers
  tuple<vector<u32>, vector<u32>, vector<u8>, vector<u32>, vector<u32>,
        vector<u8>>
      polygons;
  bool polygons_calculated = false;
  auto calculate_polygon =
      [&element_data, &polygons, &calculated =
       polygons_calculated ]() -> const auto & {
    if (!calculated) {
      auto &[nodes, element_sizes, element_indices, node_numbers,
             element_numbers] = element_data;

      polygons = get_polygon<true>(element_sizes, element_indices, node_numbers,
                                   element_numbers);

      calculated = true;
    }
    return polygons;
  };

  // indices_of_nodes_on_boundary
  set<u32> nodes_on_boundary;
  bool nodes_on_boundary_calculated = false;
  auto calculate_node_on_boundary = [
    &patches, &element_data, &nodes_on_boundary, &calculate_polygon,
    &calculated = nodes_on_boundary_calculated
  ]() -> const auto & {
    if (!calculated) {
      auto &[nodes, element_sizes, element_indices, node_numbers_map,
             element_numbers_map] = element_data;
      auto &[polygon_sizes, polygon_indices, _, node_element_numbers,
             polygon_element_numbers, surface_numbers] = calculate_polygon();
      nodes_on_boundary = node_on_boundary(patches, nodes);
      calculated = true;
    }
    return nodes_on_boundary;
  };

  // sizes_of_polygon_vertices_on_boundary,
  // indices_of_polygon_vertices_on_boundary,
  // numbers_of_nodes_of_polygon_on_boundary
  // numbers_of_element_whose_polygon_on_boundary
  // numbers_of_surface_on_boundary
  tuple<vector<u32>, vector<u32>, vector<u32>, vector<u32>, vector<u8>>
      polygons_on_boundary;
  bool polygons_on_boundary_calculated = false;
  auto calculate_polygon_on_boundary = [
    &element_data, &patches, &polygons_on_boundary, &calculate_polygon,
    &calculated = polygons_on_boundary_calculated
  ]() -> const auto & {
    if (!calculated) {
      auto &[nodes, element_sizes, element_indices, node_numbers,
             element_numbers] = element_data;
      auto &[polygon_sizes, polygon_indices, polygon_element_sizes,
             polygon_node_numbers, polygon_element_numbers, surface_numbers] =
          calculate_polygon();
      polygons_on_boundary = polygon_on_boundary(
          patches, nodes, polygon_sizes, polygon_indices, polygon_node_numbers,
          polygon_element_numbers, surface_numbers);
      calculated = true;
    }
    return polygons_on_boundary;
  };

  // Polygon surface number, polygon surface centroid positions and polygon
  // average
  tuple<vector<u8>, vector<float>, vector<float>> polygons_average;
  bool average_calculated = false;
  auto calculate_polygon_average = [
        &calculated = average_calculated, &polygons_average, &calculate_polygon,
        &calculate_polygon_on_boundary, &patches, &element_data, &frames
  ]() -> const auto & {
    auto &[polygon_sizes, polygon_indices, element_sizes, polygon_node_numbers,
           polygon_element_numbers, surface_numbers] = calculate_polygon();

    auto &[sizes_of_polygon_vertices_on_boundary,
           indices_of_polygon_vertices_on_boundary,
           indices_of_node_whose_polygon_on_boundary,
           indices_of_element_whose_polygon_on_boundary,
           numbers_of_surfaces_on_boundary] = calculate_polygon_on_boundary();

    if (!calculated) {
      auto &[nodes, element_sizes, element_indices, node_numbers,
             element_numbers] = element_data;
      polygons_average =
          polygon_average(patches, nodes, numbers_of_surfaces_on_boundary,
                          sizes_of_polygon_vertices_on_boundary,
                          indices_of_polygon_vertices_on_boundary, frames);
    }
    return polygons_average;
  };

  while (true) {

    if (!patches.empty())
      cout << R"(
s - Select a patch.
S - Settings.
p - Visulize patches.
P - Visulize patches and select some.
l - Select patches.
M - Visualize merged patches.)";

#if GRAPHICS_ENABLED
    if (!patches.empty() && elem_avail())
      cout << R"(
y - Visualize polygons on current patch.
n - Visualize nodes on current patch.
r - Visualize elements primitives on current patch.
Y - Visualize polygons on boundary.
N - Visualize nodes on boundary.
R - Visualize elements primitives on boundary.)";
#endif // GRAPHICS_ENABLED

    if (!patches.empty() && !frames.empty() && elem_avail())
      cout << R"(
l - Write APDL output.
A - Review output by visualizing them.)";

    cout <<
        R"(
d - Discard.
)";
    char opt = 'd';
    cin >> opt;
    switch (opt) {
    case 'd':
      return;
    case 's':
      select_patch(patches);
      break;
    case 'p':
      visualize_patch(patches);
      break;
    case 'P':
      visualize_patch_selection(patches);
      break;
    case 'S':
      mesh_settings();
      break;
#if GRAPHICS_ENABLED
    case 'y':
      if (selected_patch < patches.size() && elem_avail()) {
        auto [polygon_sizes, polygon_indices, _0, _1, _2, surface_numbers] =
            get_polygon<false>(element_sizes, element_indices, {}, {});
        auto [sizes_of_polygon_vertices_on_boundary,
              indices_of_polygon_vertices_on_boundary, _3, _4, _5] =
            polygon_on_boundary({patches[selected_patch]}, nodes, polygon_sizes,
                                polygon_indices, {}, {}, surface_numbers);

        visualize_polygons(nodes, sizes_of_polygon_vertices_on_boundary,
                           indices_of_polygon_vertices_on_boundary);
      }
      break;
    case 'n':
      if (selected_patch < patches.size() && elem_avail()) {
        auto indices_of_node_on_boundary =
            node_on_boundary({patches[selected_patch]}, nodes);

        visualize_nodes(nodes, vector<u32>(indices_of_node_on_boundary.begin(),
                                           indices_of_node_on_boundary.end()));
      }
      break;
    case 'r':
      if (selected_patch < patches.size() && elem_avail()) {
        auto primitive_indices =
            get_primitives(nodes, element_sizes, element_indices, wireframe);
        auto primitive_indices_on_boundary = primitive_on_boundary(
            {patches[selected_patch]}, nodes, primitive_indices, wireframe);
        visualize_primitives(nodes, primitive_indices_on_boundary);
      }
      break;
    case 'Y':
      if (!patches.empty() && elem_avail()) {
        auto &[polygon_sizes, polygon_indices, _0, _1, _2, _3] =
            calculate_polygon();
        auto &[sizes_of_polygon_vertices_on_boundary,
               indices_of_polygon_vertices_on_boundary,
               indices_of_node_whose_polygon_on_boundary,
               indices_of_element_whose_polygon_on_boundary, surface_numbers] =
            calculate_polygon_on_boundary();

        visualize_polygons(nodes, sizes_of_polygon_vertices_on_boundary,
                           indices_of_polygon_vertices_on_boundary);
      }
      break;
    case 'N':
      if (!patches.empty() && elem_avail()) {
        auto indices_of_node_on_boundary = node_on_boundary(patches, nodes);

        visualize_nodes(nodes, vector<u32>(indices_of_node_on_boundary.begin(),
                                           indices_of_node_on_boundary.end()));
      }
      break;
    case 'R':
      if (!patches.empty() && elem_avail()) {
        auto primitive_indices =
            get_primitives(nodes, element_sizes, element_indices, wireframe);
        auto primitive_indices_on_boundary =
            primitive_on_boundary(patches, nodes, primitive_indices, wireframe);
        visualize_primitives(nodes, primitive_indices_on_boundary);
      }
      break;
    case 'A':
      if (!patches.empty() && elem_avail()) {
        auto &[_, centroid, boundary_data] = calculate_polygon_average();
        if (centroid.empty() || boundary_data.empty()) {
          cerr << "No surfaces on boundary found\n";
          break;
        }
        assert(frames.size() <= numeric_limits<u32>::max());
        visualize_nodes(centroid, boundary_data, (u32)frames.size());
      }
      break;
    case 'M':
      if (!patches.empty()) {
        visualize_regions(merge(patches));
      }
      break;
#endif // GRAPHICS_ENABLED
    case 'l':
      if (!patches.empty() && !frames.empty() && elem_avail()) {
        auto opt1 = request_file_by_name([](const path &p) { return true; },
                                         "APDL output");
        auto opt2 = request_file_by_name([](const path &p) { return true; },
                                         "APDL output directory");
        if (!opt1)
          break;
        auto out = ofstream(opt1.value());
        if (!out) {
          FILE_OPEN_FAILED;
        }
        if (!opt2)
          break;
        auto &dir = opt2.value();
        if (!exists(dir)) {
          auto suc = create_directory(dir);
          if (!suc) {
            DIRECTORY_CREATE_FAILED;
          }
        }
        out << "/PREP7" << endl;

        const auto N = frames.size();

        auto &[on_boundary_vertex_sizes, on_boundary_vertex_indices,
               node_numbers, element_numbers, surface_numbers] =
            calculate_polygon_on_boundary();
        auto &[on_boundary_surface_numbers, _2, boundary_data] =
            calculate_polygon_average();

        assert(on_boundary_surface_numbers.size() ==
               on_boundary_vertex_sizes.size());

        vector<vector<u32>::const_iterator> ps;
        ps.reserve(on_boundary_surface_numbers.size() + 1);
        vector<vector<float>::const_iterator> Ps;
        Ps.reserve(on_boundary_surface_numbers.size() + 1);
        {
          auto p = on_boundary_vertex_indices.begin();
          auto e = on_boundary_vertex_indices.end();
          auto P = boundary_data.begin();
          auto E = boundary_data.end();

          for (auto sz : on_boundary_vertex_sizes) {
            assert(p < e);
            assert(P < E);
            ps.push_back(p);
            Ps.push_back(P);
            p += sz;
            P += N;
          }
          ps.push_back(p);
          Ps.push_back(P);
          assert(p == e);
          assert(P == E);
        }

        assert(element_numbers.size() == on_boundary_surface_numbers.size());
#pragma omp parallel for schedule(dynamic, 5000)
        for (long i = 0; i < on_boundary_vertex_sizes.size(); ++i) {
          write_table(out, opt2.value(), "HFLUX", element_numbers[i],
                      on_boundary_surface_numbers[i], ps[i], ps[i + 1], frames,
                      Ps[i], Ps[i + 1]);
        }
        out << "FINISH" << endl
            << "/SOLU" << endl
            << "ALLSEL,ALL" << endl
            << "/NERR,,99999999" << endl;
      }
      break;
    default:
      COMMAND_NOT_FOUND;
      break;
    }
  }
}