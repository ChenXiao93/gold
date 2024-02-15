import os
from moviepy.editor import *
from pydub import AudioSegment

AudioSegment.converter = r"D:\\AppGallery\\Software\\ffmpeg\\ffmpeg-5.0.1-essentials_build\\bin\\ffmpeg.exe"

# def match_target_amplitude(sound, target_dBFS):
#     change_in_dBFS = target_dBFS - sound.dBFS
#     return( sound.apply_gain(change_in_dBFS))
# audio = AudioSegment.from_mp3('./.mp3')


def Download():

    link = ["https://www.bilibili.com/video/BV1Z8411v7e5"]

    for i in link:
        # os.system('you-get -o D:\Python %s' % i) # 下载视频
        # os.system('you-get -i %s' % i)  # 查看当前视频的清晰度和格式
        os.system('you-get --format=dash-flv480 -o D:\视频下载\ %s' % i)  # 下载视频

    return ()


def convert():
    link = ["二次火爆全网的《桃花诺》，小姐姐开口就惊了！"]
    for i in range(len(link)):
        video = VideoFileClip("D:\视频下载\%s.mp4" % link[i])
        audio = video.audio
        # audio.write_audiofile('test.wav')
        audio.write_audiofile("D:\视频下载\%s.mp3" % link[i])
    return ()

##################################################
Download()
# convert()


# song = AudioSegment.from_mp3(r"D:\\Python\\est.mp3")










