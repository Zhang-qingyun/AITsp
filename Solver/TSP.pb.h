// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: TSP.proto

#ifndef PROTOBUF_INCLUDED_TSP_2eproto
#define PROTOBUF_INCLUDED_TSP_2eproto

#include <string>

#include <google/protobuf/stubs/common.h>

#if GOOGLE_PROTOBUF_VERSION < 3006001
#error This file was generated by a newer version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please update
#error your headers.
#endif
#if 3006001 < GOOGLE_PROTOBUF_MIN_PROTOC_VERSION
#error This file was generated by an older version of protoc which is
#error incompatible with your Protocol Buffer headers.  Please
#error regenerate this file with a newer version of protoc.
#endif

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/arena.h>
#include <google/protobuf/arenastring.h>
#include <google/protobuf/generated_message_table_driven.h>
#include <google/protobuf/generated_message_util.h>
#include <google/protobuf/inlined_string_field.h>
#include <google/protobuf/metadata.h>
#include <google/protobuf/message.h>
#include <google/protobuf/repeated_field.h>  // IWYU pragma: export
#include <google/protobuf/extension_set.h>  // IWYU pragma: export
#include <google/protobuf/unknown_field_set.h>
// @@protoc_insertion_point(includes)
#define PROTOBUF_INTERNAL_EXPORT_protobuf_TSP_2eproto 

namespace protobuf_TSP_2eproto {
// Internal implementation detail -- do not use these members.
struct TableStruct {
  static const ::google::protobuf::internal::ParseTableField entries[];
  static const ::google::protobuf::internal::AuxillaryParseTableField aux[];
  static const ::google::protobuf::internal::ParseTable schema[6];
  static const ::google::protobuf::internal::FieldMetadata field_metadata[];
  static const ::google::protobuf::internal::SerializationTable serialization_table[];
  static const ::google::protobuf::uint32 offsets[];
};
void AddDescriptors();
}  // namespace protobuf_TSP_2eproto
namespace pb {
class TSP;
class TSPDefaultTypeInternal;
extern TSPDefaultTypeInternal _TSP_default_instance_;
class TSP_Edge;
class TSP_EdgeDefaultTypeInternal;
extern TSP_EdgeDefaultTypeInternal _TSP_Edge_default_instance_;
class TSP_Input;
class TSP_InputDefaultTypeInternal;
extern TSP_InputDefaultTypeInternal _TSP_Input_default_instance_;
class TSP_Node;
class TSP_NodeDefaultTypeInternal;
extern TSP_NodeDefaultTypeInternal _TSP_Node_default_instance_;
class TSP_Output;
class TSP_OutputDefaultTypeInternal;
extern TSP_OutputDefaultTypeInternal _TSP_Output_default_instance_;
class TSP_UndirectGraph;
class TSP_UndirectGraphDefaultTypeInternal;
extern TSP_UndirectGraphDefaultTypeInternal _TSP_UndirectGraph_default_instance_;
}  // namespace pb
namespace google {
namespace protobuf {
template<> ::pb::TSP* Arena::CreateMaybeMessage<::pb::TSP>(Arena*);
template<> ::pb::TSP_Edge* Arena::CreateMaybeMessage<::pb::TSP_Edge>(Arena*);
template<> ::pb::TSP_Input* Arena::CreateMaybeMessage<::pb::TSP_Input>(Arena*);
template<> ::pb::TSP_Node* Arena::CreateMaybeMessage<::pb::TSP_Node>(Arena*);
template<> ::pb::TSP_Output* Arena::CreateMaybeMessage<::pb::TSP_Output>(Arena*);
template<> ::pb::TSP_UndirectGraph* Arena::CreateMaybeMessage<::pb::TSP_UndirectGraph>(Arena*);
}  // namespace protobuf
}  // namespace google
namespace pb {

// ===================================================================

class TSP_Input : public ::google::protobuf::Message /* @@protoc_insertion_point(class_definition:pb.TSP.Input) */ {
 public:
  TSP_Input();
  virtual ~TSP_Input();

  TSP_Input(const TSP_Input& from);

  inline TSP_Input& operator=(const TSP_Input& from) {
    CopyFrom(from);
    return *this;
  }
  #if LANG_CXX11
  TSP_Input(TSP_Input&& from) noexcept
    : TSP_Input() {
    *this = ::std::move(from);
  }

  inline TSP_Input& operator=(TSP_Input&& from) noexcept {
    if (GetArenaNoVirtual() == from.GetArenaNoVirtual()) {
      if (this != &from) InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }
  #endif
  static const ::google::protobuf::Descriptor* descriptor();
  static const TSP_Input& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const TSP_Input* internal_default_instance() {
    return reinterpret_cast<const TSP_Input*>(
               &_TSP_Input_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  void Swap(TSP_Input* other);
  friend void swap(TSP_Input& a, TSP_Input& b) {
    a.Swap(&b);
  }

  // implements Message ----------------------------------------------

  inline TSP_Input* New() const final {
    return CreateMaybeMessage<TSP_Input>(NULL);
  }

  TSP_Input* New(::google::protobuf::Arena* arena) const final {
    return CreateMaybeMessage<TSP_Input>(arena);
  }
  void CopyFrom(const ::google::protobuf::Message& from) final;
  void MergeFrom(const ::google::protobuf::Message& from) final;
  void CopyFrom(const TSP_Input& from);
  void MergeFrom(const TSP_Input& from);
  void Clear() final;
  bool IsInitialized() const final;

  size_t ByteSizeLong() const final;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input) final;
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const final;
  ::google::protobuf::uint8* InternalSerializeWithCachedSizesToArray(
      bool deterministic, ::google::protobuf::uint8* target) const final;
  int GetCachedSize() const final { return _cached_size_.Get(); }

  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const final;
  void InternalSwap(TSP_Input* other);
  private:
  inline ::google::protobuf::Arena* GetArenaNoVirtual() const {
    return NULL;
  }
  inline void* MaybeArenaPtr() const {
    return NULL;
  }
  public:

  ::google::protobuf::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  // .pb.TSP.UndirectGraph graph = 1;
  bool has_graph() const;
  void clear_graph();
  static const int kGraphFieldNumber = 1;
  private:
  const ::pb::TSP_UndirectGraph& _internal_graph() const;
  public:
  const ::pb::TSP_UndirectGraph& graph() const;
  ::pb::TSP_UndirectGraph* release_graph();
  ::pb::TSP_UndirectGraph* mutable_graph();
  void set_allocated_graph(::pb::TSP_UndirectGraph* graph);

  // int32 centerNum = 2;
  void clear_centernum();
  static const int kCenterNumFieldNumber = 2;
  ::google::protobuf::int32 centernum() const;
  void set_centernum(::google::protobuf::int32 value);

  // @@protoc_insertion_point(class_scope:pb.TSP.Input)
 private:

  ::google::protobuf::internal::InternalMetadataWithArena _internal_metadata_;
  ::pb::TSP_UndirectGraph* graph_;
  ::google::protobuf::int32 centernum_;
  mutable ::google::protobuf::internal::CachedSize _cached_size_;
  friend struct ::protobuf_TSP_2eproto::TableStruct;
};
// -------------------------------------------------------------------

class TSP_Output : public ::google::protobuf::Message /* @@protoc_insertion_point(class_definition:pb.TSP.Output) */ {
 public:
  TSP_Output();
  virtual ~TSP_Output();

  TSP_Output(const TSP_Output& from);

  inline TSP_Output& operator=(const TSP_Output& from) {
    CopyFrom(from);
    return *this;
  }
  #if LANG_CXX11
  TSP_Output(TSP_Output&& from) noexcept
    : TSP_Output() {
    *this = ::std::move(from);
  }

  inline TSP_Output& operator=(TSP_Output&& from) noexcept {
    if (GetArenaNoVirtual() == from.GetArenaNoVirtual()) {
      if (this != &from) InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }
  #endif
  static const ::google::protobuf::Descriptor* descriptor();
  static const TSP_Output& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const TSP_Output* internal_default_instance() {
    return reinterpret_cast<const TSP_Output*>(
               &_TSP_Output_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    1;

  void Swap(TSP_Output* other);
  friend void swap(TSP_Output& a, TSP_Output& b) {
    a.Swap(&b);
  }

  // implements Message ----------------------------------------------

  inline TSP_Output* New() const final {
    return CreateMaybeMessage<TSP_Output>(NULL);
  }

  TSP_Output* New(::google::protobuf::Arena* arena) const final {
    return CreateMaybeMessage<TSP_Output>(arena);
  }
  void CopyFrom(const ::google::protobuf::Message& from) final;
  void MergeFrom(const ::google::protobuf::Message& from) final;
  void CopyFrom(const TSP_Output& from);
  void MergeFrom(const TSP_Output& from);
  void Clear() final;
  bool IsInitialized() const final;

  size_t ByteSizeLong() const final;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input) final;
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const final;
  ::google::protobuf::uint8* InternalSerializeWithCachedSizesToArray(
      bool deterministic, ::google::protobuf::uint8* target) const final;
  int GetCachedSize() const final { return _cached_size_.Get(); }

  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const final;
  void InternalSwap(TSP_Output* other);
  private:
  inline ::google::protobuf::Arena* GetArenaNoVirtual() const {
    return NULL;
  }
  inline void* MaybeArenaPtr() const {
    return NULL;
  }
  public:

  ::google::protobuf::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  // repeated int32 paths = 1;
  int paths_size() const;
  void clear_paths();
  static const int kPathsFieldNumber = 1;
  ::google::protobuf::int32 paths(int index) const;
  void set_paths(int index, ::google::protobuf::int32 value);
  void add_paths(::google::protobuf::int32 value);
  const ::google::protobuf::RepeatedField< ::google::protobuf::int32 >&
      paths() const;
  ::google::protobuf::RepeatedField< ::google::protobuf::int32 >*
      mutable_paths();

  // @@protoc_insertion_point(class_scope:pb.TSP.Output)
 private:

  ::google::protobuf::internal::InternalMetadataWithArena _internal_metadata_;
  ::google::protobuf::RepeatedField< ::google::protobuf::int32 > paths_;
  mutable int _paths_cached_byte_size_;
  mutable ::google::protobuf::internal::CachedSize _cached_size_;
  friend struct ::protobuf_TSP_2eproto::TableStruct;
};
// -------------------------------------------------------------------

class TSP_UndirectGraph : public ::google::protobuf::Message /* @@protoc_insertion_point(class_definition:pb.TSP.UndirectGraph) */ {
 public:
  TSP_UndirectGraph();
  virtual ~TSP_UndirectGraph();

  TSP_UndirectGraph(const TSP_UndirectGraph& from);

  inline TSP_UndirectGraph& operator=(const TSP_UndirectGraph& from) {
    CopyFrom(from);
    return *this;
  }
  #if LANG_CXX11
  TSP_UndirectGraph(TSP_UndirectGraph&& from) noexcept
    : TSP_UndirectGraph() {
    *this = ::std::move(from);
  }

  inline TSP_UndirectGraph& operator=(TSP_UndirectGraph&& from) noexcept {
    if (GetArenaNoVirtual() == from.GetArenaNoVirtual()) {
      if (this != &from) InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }
  #endif
  static const ::google::protobuf::Descriptor* descriptor();
  static const TSP_UndirectGraph& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const TSP_UndirectGraph* internal_default_instance() {
    return reinterpret_cast<const TSP_UndirectGraph*>(
               &_TSP_UndirectGraph_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    2;

  void Swap(TSP_UndirectGraph* other);
  friend void swap(TSP_UndirectGraph& a, TSP_UndirectGraph& b) {
    a.Swap(&b);
  }

  // implements Message ----------------------------------------------

  inline TSP_UndirectGraph* New() const final {
    return CreateMaybeMessage<TSP_UndirectGraph>(NULL);
  }

  TSP_UndirectGraph* New(::google::protobuf::Arena* arena) const final {
    return CreateMaybeMessage<TSP_UndirectGraph>(arena);
  }
  void CopyFrom(const ::google::protobuf::Message& from) final;
  void MergeFrom(const ::google::protobuf::Message& from) final;
  void CopyFrom(const TSP_UndirectGraph& from);
  void MergeFrom(const TSP_UndirectGraph& from);
  void Clear() final;
  bool IsInitialized() const final;

  size_t ByteSizeLong() const final;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input) final;
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const final;
  ::google::protobuf::uint8* InternalSerializeWithCachedSizesToArray(
      bool deterministic, ::google::protobuf::uint8* target) const final;
  int GetCachedSize() const final { return _cached_size_.Get(); }

  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const final;
  void InternalSwap(TSP_UndirectGraph* other);
  private:
  inline ::google::protobuf::Arena* GetArenaNoVirtual() const {
    return NULL;
  }
  inline void* MaybeArenaPtr() const {
    return NULL;
  }
  public:

  ::google::protobuf::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  // repeated .pb.TSP.Edge edges = 1;
  int edges_size() const;
  void clear_edges();
  static const int kEdgesFieldNumber = 1;
  ::pb::TSP_Edge* mutable_edges(int index);
  ::google::protobuf::RepeatedPtrField< ::pb::TSP_Edge >*
      mutable_edges();
  const ::pb::TSP_Edge& edges(int index) const;
  ::pb::TSP_Edge* add_edges();
  const ::google::protobuf::RepeatedPtrField< ::pb::TSP_Edge >&
      edges() const;

  // repeated .pb.TSP.Node nodes = 3;
  int nodes_size() const;
  void clear_nodes();
  static const int kNodesFieldNumber = 3;
  ::pb::TSP_Node* mutable_nodes(int index);
  ::google::protobuf::RepeatedPtrField< ::pb::TSP_Node >*
      mutable_nodes();
  const ::pb::TSP_Node& nodes(int index) const;
  ::pb::TSP_Node* add_nodes();
  const ::google::protobuf::RepeatedPtrField< ::pb::TSP_Node >&
      nodes() const;

  // int32 nodeNum = 2;
  void clear_nodenum();
  static const int kNodeNumFieldNumber = 2;
  ::google::protobuf::int32 nodenum() const;
  void set_nodenum(::google::protobuf::int32 value);

  // @@protoc_insertion_point(class_scope:pb.TSP.UndirectGraph)
 private:

  ::google::protobuf::internal::InternalMetadataWithArena _internal_metadata_;
  ::google::protobuf::RepeatedPtrField< ::pb::TSP_Edge > edges_;
  ::google::protobuf::RepeatedPtrField< ::pb::TSP_Node > nodes_;
  ::google::protobuf::int32 nodenum_;
  mutable ::google::protobuf::internal::CachedSize _cached_size_;
  friend struct ::protobuf_TSP_2eproto::TableStruct;
};
// -------------------------------------------------------------------

class TSP_Edge : public ::google::protobuf::Message /* @@protoc_insertion_point(class_definition:pb.TSP.Edge) */ {
 public:
  TSP_Edge();
  virtual ~TSP_Edge();

  TSP_Edge(const TSP_Edge& from);

  inline TSP_Edge& operator=(const TSP_Edge& from) {
    CopyFrom(from);
    return *this;
  }
  #if LANG_CXX11
  TSP_Edge(TSP_Edge&& from) noexcept
    : TSP_Edge() {
    *this = ::std::move(from);
  }

  inline TSP_Edge& operator=(TSP_Edge&& from) noexcept {
    if (GetArenaNoVirtual() == from.GetArenaNoVirtual()) {
      if (this != &from) InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }
  #endif
  static const ::google::protobuf::Descriptor* descriptor();
  static const TSP_Edge& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const TSP_Edge* internal_default_instance() {
    return reinterpret_cast<const TSP_Edge*>(
               &_TSP_Edge_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    3;

  void Swap(TSP_Edge* other);
  friend void swap(TSP_Edge& a, TSP_Edge& b) {
    a.Swap(&b);
  }

  // implements Message ----------------------------------------------

  inline TSP_Edge* New() const final {
    return CreateMaybeMessage<TSP_Edge>(NULL);
  }

  TSP_Edge* New(::google::protobuf::Arena* arena) const final {
    return CreateMaybeMessage<TSP_Edge>(arena);
  }
  void CopyFrom(const ::google::protobuf::Message& from) final;
  void MergeFrom(const ::google::protobuf::Message& from) final;
  void CopyFrom(const TSP_Edge& from);
  void MergeFrom(const TSP_Edge& from);
  void Clear() final;
  bool IsInitialized() const final;

  size_t ByteSizeLong() const final;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input) final;
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const final;
  ::google::protobuf::uint8* InternalSerializeWithCachedSizesToArray(
      bool deterministic, ::google::protobuf::uint8* target) const final;
  int GetCachedSize() const final { return _cached_size_.Get(); }

  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const final;
  void InternalSwap(TSP_Edge* other);
  private:
  inline ::google::protobuf::Arena* GetArenaNoVirtual() const {
    return NULL;
  }
  inline void* MaybeArenaPtr() const {
    return NULL;
  }
  public:

  ::google::protobuf::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  // int32 source = 1;
  void clear_source();
  static const int kSourceFieldNumber = 1;
  ::google::protobuf::int32 source() const;
  void set_source(::google::protobuf::int32 value);

  // int32 target = 2;
  void clear_target();
  static const int kTargetFieldNumber = 2;
  ::google::protobuf::int32 target() const;
  void set_target(::google::protobuf::int32 value);

  // int32 length = 3;
  void clear_length();
  static const int kLengthFieldNumber = 3;
  ::google::protobuf::int32 length() const;
  void set_length(::google::protobuf::int32 value);

  // @@protoc_insertion_point(class_scope:pb.TSP.Edge)
 private:

  ::google::protobuf::internal::InternalMetadataWithArena _internal_metadata_;
  ::google::protobuf::int32 source_;
  ::google::protobuf::int32 target_;
  ::google::protobuf::int32 length_;
  mutable ::google::protobuf::internal::CachedSize _cached_size_;
  friend struct ::protobuf_TSP_2eproto::TableStruct;
};
// -------------------------------------------------------------------

class TSP_Node : public ::google::protobuf::Message /* @@protoc_insertion_point(class_definition:pb.TSP.Node) */ {
 public:
  TSP_Node();
  virtual ~TSP_Node();

  TSP_Node(const TSP_Node& from);

  inline TSP_Node& operator=(const TSP_Node& from) {
    CopyFrom(from);
    return *this;
  }
  #if LANG_CXX11
  TSP_Node(TSP_Node&& from) noexcept
    : TSP_Node() {
    *this = ::std::move(from);
  }

  inline TSP_Node& operator=(TSP_Node&& from) noexcept {
    if (GetArenaNoVirtual() == from.GetArenaNoVirtual()) {
      if (this != &from) InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }
  #endif
  static const ::google::protobuf::Descriptor* descriptor();
  static const TSP_Node& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const TSP_Node* internal_default_instance() {
    return reinterpret_cast<const TSP_Node*>(
               &_TSP_Node_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    4;

  void Swap(TSP_Node* other);
  friend void swap(TSP_Node& a, TSP_Node& b) {
    a.Swap(&b);
  }

  // implements Message ----------------------------------------------

  inline TSP_Node* New() const final {
    return CreateMaybeMessage<TSP_Node>(NULL);
  }

  TSP_Node* New(::google::protobuf::Arena* arena) const final {
    return CreateMaybeMessage<TSP_Node>(arena);
  }
  void CopyFrom(const ::google::protobuf::Message& from) final;
  void MergeFrom(const ::google::protobuf::Message& from) final;
  void CopyFrom(const TSP_Node& from);
  void MergeFrom(const TSP_Node& from);
  void Clear() final;
  bool IsInitialized() const final;

  size_t ByteSizeLong() const final;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input) final;
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const final;
  ::google::protobuf::uint8* InternalSerializeWithCachedSizesToArray(
      bool deterministic, ::google::protobuf::uint8* target) const final;
  int GetCachedSize() const final { return _cached_size_.Get(); }

  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const final;
  void InternalSwap(TSP_Node* other);
  private:
  inline ::google::protobuf::Arena* GetArenaNoVirtual() const {
    return NULL;
  }
  inline void* MaybeArenaPtr() const {
    return NULL;
  }
  public:

  ::google::protobuf::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  // double x = 1;
  void clear_x();
  static const int kXFieldNumber = 1;
  double x() const;
  void set_x(double value);

  // double y = 2;
  void clear_y();
  static const int kYFieldNumber = 2;
  double y() const;
  void set_y(double value);

  // @@protoc_insertion_point(class_scope:pb.TSP.Node)
 private:

  ::google::protobuf::internal::InternalMetadataWithArena _internal_metadata_;
  double x_;
  double y_;
  mutable ::google::protobuf::internal::CachedSize _cached_size_;
  friend struct ::protobuf_TSP_2eproto::TableStruct;
};
// -------------------------------------------------------------------

class TSP : public ::google::protobuf::Message /* @@protoc_insertion_point(class_definition:pb.TSP) */ {
 public:
  TSP();
  virtual ~TSP();

  TSP(const TSP& from);

  inline TSP& operator=(const TSP& from) {
    CopyFrom(from);
    return *this;
  }
  #if LANG_CXX11
  TSP(TSP&& from) noexcept
    : TSP() {
    *this = ::std::move(from);
  }

  inline TSP& operator=(TSP&& from) noexcept {
    if (GetArenaNoVirtual() == from.GetArenaNoVirtual()) {
      if (this != &from) InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }
  #endif
  static const ::google::protobuf::Descriptor* descriptor();
  static const TSP& default_instance();

  static void InitAsDefaultInstance();  // FOR INTERNAL USE ONLY
  static inline const TSP* internal_default_instance() {
    return reinterpret_cast<const TSP*>(
               &_TSP_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    5;

  void Swap(TSP* other);
  friend void swap(TSP& a, TSP& b) {
    a.Swap(&b);
  }

  // implements Message ----------------------------------------------

  inline TSP* New() const final {
    return CreateMaybeMessage<TSP>(NULL);
  }

  TSP* New(::google::protobuf::Arena* arena) const final {
    return CreateMaybeMessage<TSP>(arena);
  }
  void CopyFrom(const ::google::protobuf::Message& from) final;
  void MergeFrom(const ::google::protobuf::Message& from) final;
  void CopyFrom(const TSP& from);
  void MergeFrom(const TSP& from);
  void Clear() final;
  bool IsInitialized() const final;

  size_t ByteSizeLong() const final;
  bool MergePartialFromCodedStream(
      ::google::protobuf::io::CodedInputStream* input) final;
  void SerializeWithCachedSizes(
      ::google::protobuf::io::CodedOutputStream* output) const final;
  ::google::protobuf::uint8* InternalSerializeWithCachedSizesToArray(
      bool deterministic, ::google::protobuf::uint8* target) const final;
  int GetCachedSize() const final { return _cached_size_.Get(); }

  private:
  void SharedCtor();
  void SharedDtor();
  void SetCachedSize(int size) const final;
  void InternalSwap(TSP* other);
  private:
  inline ::google::protobuf::Arena* GetArenaNoVirtual() const {
    return NULL;
  }
  inline void* MaybeArenaPtr() const {
    return NULL;
  }
  public:

  ::google::protobuf::Metadata GetMetadata() const final;

  // nested types ----------------------------------------------------

  typedef TSP_Input Input;
  typedef TSP_Output Output;
  typedef TSP_UndirectGraph UndirectGraph;
  typedef TSP_Edge Edge;
  typedef TSP_Node Node;

  // accessors -------------------------------------------------------

  // @@protoc_insertion_point(class_scope:pb.TSP)
 private:

  ::google::protobuf::internal::InternalMetadataWithArena _internal_metadata_;
  mutable ::google::protobuf::internal::CachedSize _cached_size_;
  friend struct ::protobuf_TSP_2eproto::TableStruct;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// TSP_Input

// .pb.TSP.UndirectGraph graph = 1;
inline bool TSP_Input::has_graph() const {
  return this != internal_default_instance() && graph_ != NULL;
}
inline void TSP_Input::clear_graph() {
  if (GetArenaNoVirtual() == NULL && graph_ != NULL) {
    delete graph_;
  }
  graph_ = NULL;
}
inline const ::pb::TSP_UndirectGraph& TSP_Input::_internal_graph() const {
  return *graph_;
}
inline const ::pb::TSP_UndirectGraph& TSP_Input::graph() const {
  const ::pb::TSP_UndirectGraph* p = graph_;
  // @@protoc_insertion_point(field_get:pb.TSP.Input.graph)
  return p != NULL ? *p : *reinterpret_cast<const ::pb::TSP_UndirectGraph*>(
      &::pb::_TSP_UndirectGraph_default_instance_);
}
inline ::pb::TSP_UndirectGraph* TSP_Input::release_graph() {
  // @@protoc_insertion_point(field_release:pb.TSP.Input.graph)
  
  ::pb::TSP_UndirectGraph* temp = graph_;
  graph_ = NULL;
  return temp;
}
inline ::pb::TSP_UndirectGraph* TSP_Input::mutable_graph() {
  
  if (graph_ == NULL) {
    auto* p = CreateMaybeMessage<::pb::TSP_UndirectGraph>(GetArenaNoVirtual());
    graph_ = p;
  }
  // @@protoc_insertion_point(field_mutable:pb.TSP.Input.graph)
  return graph_;
}
inline void TSP_Input::set_allocated_graph(::pb::TSP_UndirectGraph* graph) {
  ::google::protobuf::Arena* message_arena = GetArenaNoVirtual();
  if (message_arena == NULL) {
    delete graph_;
  }
  if (graph) {
    ::google::protobuf::Arena* submessage_arena = NULL;
    if (message_arena != submessage_arena) {
      graph = ::google::protobuf::internal::GetOwnedMessage(
          message_arena, graph, submessage_arena);
    }
    
  } else {
    
  }
  graph_ = graph;
  // @@protoc_insertion_point(field_set_allocated:pb.TSP.Input.graph)
}

// int32 centerNum = 2;
inline void TSP_Input::clear_centernum() {
  centernum_ = 0;
}
inline ::google::protobuf::int32 TSP_Input::centernum() const {
  // @@protoc_insertion_point(field_get:pb.TSP.Input.centerNum)
  return centernum_;
}
inline void TSP_Input::set_centernum(::google::protobuf::int32 value) {
  
  centernum_ = value;
  // @@protoc_insertion_point(field_set:pb.TSP.Input.centerNum)
}

// -------------------------------------------------------------------

// TSP_Output

// repeated int32 paths = 1;
inline int TSP_Output::paths_size() const {
  return paths_.size();
}
inline void TSP_Output::clear_paths() {
  paths_.Clear();
}
inline ::google::protobuf::int32 TSP_Output::paths(int index) const {
  // @@protoc_insertion_point(field_get:pb.TSP.Output.paths)
  return paths_.Get(index);
}
inline void TSP_Output::set_paths(int index, ::google::protobuf::int32 value) {
  paths_.Set(index, value);
  // @@protoc_insertion_point(field_set:pb.TSP.Output.paths)
}
inline void TSP_Output::add_paths(::google::protobuf::int32 value) {
  paths_.Add(value);
  // @@protoc_insertion_point(field_add:pb.TSP.Output.paths)
}
inline const ::google::protobuf::RepeatedField< ::google::protobuf::int32 >&
TSP_Output::paths() const {
  // @@protoc_insertion_point(field_list:pb.TSP.Output.paths)
  return paths_;
}
inline ::google::protobuf::RepeatedField< ::google::protobuf::int32 >*
TSP_Output::mutable_paths() {
  // @@protoc_insertion_point(field_mutable_list:pb.TSP.Output.paths)
  return &paths_;
}

// -------------------------------------------------------------------

// TSP_UndirectGraph

// repeated .pb.TSP.Edge edges = 1;
inline int TSP_UndirectGraph::edges_size() const {
  return edges_.size();
}
inline void TSP_UndirectGraph::clear_edges() {
  edges_.Clear();
}
inline ::pb::TSP_Edge* TSP_UndirectGraph::mutable_edges(int index) {
  // @@protoc_insertion_point(field_mutable:pb.TSP.UndirectGraph.edges)
  return edges_.Mutable(index);
}
inline ::google::protobuf::RepeatedPtrField< ::pb::TSP_Edge >*
TSP_UndirectGraph::mutable_edges() {
  // @@protoc_insertion_point(field_mutable_list:pb.TSP.UndirectGraph.edges)
  return &edges_;
}
inline const ::pb::TSP_Edge& TSP_UndirectGraph::edges(int index) const {
  // @@protoc_insertion_point(field_get:pb.TSP.UndirectGraph.edges)
  return edges_.Get(index);
}
inline ::pb::TSP_Edge* TSP_UndirectGraph::add_edges() {
  // @@protoc_insertion_point(field_add:pb.TSP.UndirectGraph.edges)
  return edges_.Add();
}
inline const ::google::protobuf::RepeatedPtrField< ::pb::TSP_Edge >&
TSP_UndirectGraph::edges() const {
  // @@protoc_insertion_point(field_list:pb.TSP.UndirectGraph.edges)
  return edges_;
}

// int32 nodeNum = 2;
inline void TSP_UndirectGraph::clear_nodenum() {
  nodenum_ = 0;
}
inline ::google::protobuf::int32 TSP_UndirectGraph::nodenum() const {
  // @@protoc_insertion_point(field_get:pb.TSP.UndirectGraph.nodeNum)
  return nodenum_;
}
inline void TSP_UndirectGraph::set_nodenum(::google::protobuf::int32 value) {
  
  nodenum_ = value;
  // @@protoc_insertion_point(field_set:pb.TSP.UndirectGraph.nodeNum)
}

// repeated .pb.TSP.Node nodes = 3;
inline int TSP_UndirectGraph::nodes_size() const {
  return nodes_.size();
}
inline void TSP_UndirectGraph::clear_nodes() {
  nodes_.Clear();
}
inline ::pb::TSP_Node* TSP_UndirectGraph::mutable_nodes(int index) {
  // @@protoc_insertion_point(field_mutable:pb.TSP.UndirectGraph.nodes)
  return nodes_.Mutable(index);
}
inline ::google::protobuf::RepeatedPtrField< ::pb::TSP_Node >*
TSP_UndirectGraph::mutable_nodes() {
  // @@protoc_insertion_point(field_mutable_list:pb.TSP.UndirectGraph.nodes)
  return &nodes_;
}
inline const ::pb::TSP_Node& TSP_UndirectGraph::nodes(int index) const {
  // @@protoc_insertion_point(field_get:pb.TSP.UndirectGraph.nodes)
  return nodes_.Get(index);
}
inline ::pb::TSP_Node* TSP_UndirectGraph::add_nodes() {
  // @@protoc_insertion_point(field_add:pb.TSP.UndirectGraph.nodes)
  return nodes_.Add();
}
inline const ::google::protobuf::RepeatedPtrField< ::pb::TSP_Node >&
TSP_UndirectGraph::nodes() const {
  // @@protoc_insertion_point(field_list:pb.TSP.UndirectGraph.nodes)
  return nodes_;
}

// -------------------------------------------------------------------

// TSP_Edge

// int32 source = 1;
inline void TSP_Edge::clear_source() {
  source_ = 0;
}
inline ::google::protobuf::int32 TSP_Edge::source() const {
  // @@protoc_insertion_point(field_get:pb.TSP.Edge.source)
  return source_;
}
inline void TSP_Edge::set_source(::google::protobuf::int32 value) {
  
  source_ = value;
  // @@protoc_insertion_point(field_set:pb.TSP.Edge.source)
}

// int32 target = 2;
inline void TSP_Edge::clear_target() {
  target_ = 0;
}
inline ::google::protobuf::int32 TSP_Edge::target() const {
  // @@protoc_insertion_point(field_get:pb.TSP.Edge.target)
  return target_;
}
inline void TSP_Edge::set_target(::google::protobuf::int32 value) {
  
  target_ = value;
  // @@protoc_insertion_point(field_set:pb.TSP.Edge.target)
}

// int32 length = 3;
inline void TSP_Edge::clear_length() {
  length_ = 0;
}
inline ::google::protobuf::int32 TSP_Edge::length() const {
  // @@protoc_insertion_point(field_get:pb.TSP.Edge.length)
  return length_;
}
inline void TSP_Edge::set_length(::google::protobuf::int32 value) {
  
  length_ = value;
  // @@protoc_insertion_point(field_set:pb.TSP.Edge.length)
}

// -------------------------------------------------------------------

// TSP_Node

// double x = 1;
inline void TSP_Node::clear_x() {
  x_ = 0;
}
inline double TSP_Node::x() const {
  // @@protoc_insertion_point(field_get:pb.TSP.Node.x)
  return x_;
}
inline void TSP_Node::set_x(double value) {
  
  x_ = value;
  // @@protoc_insertion_point(field_set:pb.TSP.Node.x)
}

// double y = 2;
inline void TSP_Node::clear_y() {
  y_ = 0;
}
inline double TSP_Node::y() const {
  // @@protoc_insertion_point(field_get:pb.TSP.Node.y)
  return y_;
}
inline void TSP_Node::set_y(double value) {
  
  y_ = value;
  // @@protoc_insertion_point(field_set:pb.TSP.Node.y)
}

// -------------------------------------------------------------------

// TSP

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__
// -------------------------------------------------------------------

// -------------------------------------------------------------------

// -------------------------------------------------------------------

// -------------------------------------------------------------------

// -------------------------------------------------------------------


// @@protoc_insertion_point(namespace_scope)

}  // namespace pb

// @@protoc_insertion_point(global_scope)

#endif  // PROTOBUF_INCLUDED_TSP_2eproto
