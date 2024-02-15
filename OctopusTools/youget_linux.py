import os
from moviepy.editor import *
from pydub import AudioSegment
AudioSegment.converter=r"/data117/chen/ffmpeg/ffmpeg/bin/ffmpeg"
def Download():
    link =["https://www.bilibili.com/video/BV1XM4y1D7XH/?spm_id_from=333.1007.tianma.1-1-1.click&vd_source=f8b68714a2849dc14fa918e11884604d"]
    for i in link:
        # os.system('you-get -o D:\Python %s' % i) # 下载视频
        os.system('you-get -i %s' % i)  # 查看当前视频的清晰度和格式
        # os.system('you-get --format=flv360 -o /data117/chen/ffmpeg/vedios/ %s' % i)  # 下载视频
    return ()
def convert():
    link = ["2020年好听又火的古风歌曲《谪仙》你还记得吗？"]
    for i in range(len(link)):
        video = VideoFileClip("/data117/chen/ffmpeg/vedios/%s.flv" % link[i])
        audio = video.audio
        # audio.write_audiofile('test.wav')
        audio.write_audiofile("/data117/chen/ffmpeg/vedios/%s.mp3" % link[i])
    return ()
##################################################
Download()
#
# convert()
